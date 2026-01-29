#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "conduit.hpp"
#include "conduit_blueprint.hpp"
#include "conduit_relay.hpp"

using namespace conduit;

// ==============================================
// SIMULATION PARAMETERS
// ==============================================
struct SimulationParams {
    int nx = 20;          // Grid points in x
    int ny = 15;          // Grid points in y
    int nz = 10;          // Grid points in z
    double dx = 0.1;      // Grid spacing
    double dy = 0.1;
    double dz = 0.1;
    double dt = 0.01;     // Time step
    double viscosity = 0.1; // Kinematic viscosity
    int total_iterations = 10;
    
    // Domain boundaries
    double x_min = 0.0;
    double x_max = nx * dx;
    double y_min = 0.0;
    double y_max = ny * dy;
    double z_min = 0.0;
    double z_max = nz * dz;
};

// ==============================================
// SIMULATION FIELD MANAGER
// ==============================================
class SimulationFields {
private:
    int nx, ny, nz;
    int n_total;
    
public:
    std::vector<double> pressure;
    std::vector<double> velocity_u;
    std::vector<double> velocity_v;
    std::vector<double> velocity_w;
    std::vector<double> temperature;
    std::vector<double> vorticity;
    
    SimulationFields(int nx_, int ny_, int nz_) 
        : nx(nx_), ny(ny_), nz(nz_), n_total(nx_ * ny_ * nz_) {
        
        pressure.resize(n_total, 0.0);
        velocity_u.resize(n_total, 0.0);
        velocity_v.resize(n_total, 0.0);
        velocity_w.resize(n_total, 0.0);
        temperature.resize(n_total, 0.0);
        vorticity.resize(n_total, 0.0);
    }
    
    int idx(int i, int j, int k) const {
        return i + j * nx + k * nx * ny;
    }
    
    void initialize_fields() {
        // Initialize with some interesting patterns
        for (int k = 0; k < nz; k++) {
            double z = k * (2.0 * M_PI / nz);
            for (int j = 0; j < ny; j++) {
                double y = j * (2.0 * M_PI / ny);
                for (int i = 0; i < nx; i++) {
                    double x = i * (2.0 * M_PI / nx);
                    int id = idx(i, j, k);
                    
                    // Taylor-Green vortex pattern
                    velocity_u[id] = sin(x) * cos(y) * cos(z);
                    velocity_v[id] = -cos(x) * sin(y) * cos(z);
                    velocity_w[id] = 0.0;
                    
                    // Pressure from Bernoulli
                    pressure[id] = 0.25 * (cos(2*x) + cos(2*y)) * cos(z);
                    
                    // Temperature gradient
                    temperature[id] = 300.0 + 50.0 * sin(0.5 * z);
                }
            }
        }
        
        compute_vorticity();
    }
    
    void compute_vorticity() {
        // Simplified vorticity computation
        for (int k = 1; k < nz-1; k++) {
            for (int j = 1; j < ny-1; j++) {
                for (int i = 1; i < nx-1; i++) {
                    int id = idx(i, j, k);
                    
                    // ∂v/∂x - ∂u/∂y (2D vorticity for simplicity)
                    double dv_dx = (velocity_v[idx(i+1, j, k)] - velocity_v[idx(i-1, j, k)]) / 2.0;
                    double du_dy = (velocity_u[idx(i, j+1, k)] - velocity_u[idx(i, j-1, k)]) / 2.0;
                    
                    vorticity[id] = dv_dx - du_dy;
                }
            }
        }
    }
};

// ==============================================
// SIMULATION ENGINE
// ==============================================
class SimulationEngine {
private:
    SimulationParams params;
    SimulationFields fields;
    int current_iteration;
    double current_time;
    
public:
    SimulationEngine(const SimulationParams& p) 
        : params(p), fields(p.nx, p.ny, p.nz), 
          current_iteration(0), current_time(0.0) {
        
        std::cout << "Initializing simulation..." << std::endl;
        std::cout << "Grid size: " << p.nx << " x " << p.ny << " x " << p.nz 
                  << " = " << (p.nx * p.ny * p.nz) << " points" << std::endl;
        
        fields.initialize_fields();
    }
    
    void run_iteration() {
        // Simplified Navier-Stokes-like update (not physically accurate, just for demo)
        
        // Create temporary arrays for updates
        std::vector<double> u_new = fields.velocity_u;
        std::vector<double> v_new = fields.velocity_v;
        std::vector<double> p_new = fields.pressure;
        std::vector<double> t_new = fields.temperature;
        
        // Apply simple diffusion and convection
        for (int k = 1; k < params.nz-1; k++) {
            for (int j = 1; j < params.ny-1; j++) {
                for (int i = 1; i < params.nx-1; i++) {
                    int id = fields.idx(i, j, k);
                    
                    // Laplacian for diffusion
                    double laplacian_u = 
                        (fields.velocity_u[fields.idx(i+1, j, k)] - 2*fields.velocity_u[id] + fields.velocity_u[fields.idx(i-1, j, k)])/(params.dx*params.dx) +
                        (fields.velocity_u[fields.idx(i, j+1, k)] - 2*fields.velocity_u[id] + fields.velocity_u[fields.idx(i, j-1, k)])/(params.dy*params.dy) +
                        (fields.velocity_u[fields.idx(i, j, k+1)] - 2*fields.velocity_u[id] + fields.velocity_u[fields.idx(i, j, k-1)])/(params.dz*params.dz);
                    
                    double laplacian_v = 
                        (fields.velocity_v[fields.idx(i+1, j, k)] - 2*fields.velocity_v[id] + fields.velocity_v[fields.idx(i-1, j, k)])/(params.dx*params.dx) +
                        (fields.velocity_v[fields.idx(i, j+1, k)] - 2*fields.velocity_v[id] + fields.velocity_v[fields.idx(i, j-1, k)])/(params.dy*params.dy) +
                        (fields.velocity_v[fields.idx(i, j, k+1)] - 2*fields.velocity_v[id] + fields.velocity_v[fields.idx(i, j, k-1)])/(params.dz*params.dz);
                    
                    // Convection terms (simple upwind)
                    double conv_u = -fields.velocity_u[id] * 
                        (fields.velocity_u[fields.idx(i+1, j, k)] - fields.velocity_u[fields.idx(i-1, j, k)])/(2*params.dx);
                    
                    double conv_v = -fields.velocity_v[id] * 
                        (fields.velocity_v[fields.idx(i, j+1, k)] - fields.velocity_v[fields.idx(i, j-1, k)])/(2*params.dy);
                    
                    // Update velocities
                    u_new[id] = fields.velocity_u[id] + params.dt * (
                        params.viscosity * laplacian_u + conv_u
                    );
                    
                    v_new[id] = fields.velocity_v[id] + params.dt * (
                        params.viscosity * laplacian_v + conv_v
                    );
                    
                    // Temperature diffusion
                    double laplacian_t = 
                        (fields.temperature[fields.idx(i+1, j, k)] - 2*fields.temperature[id] + fields.temperature[fields.idx(i-1, j, k)])/(params.dx*params.dx) +
                        (fields.temperature[fields.idx(i, j+1, k)] - 2*fields.temperature[id] + fields.temperature[fields.idx(i, j-1, k)])/(params.dy*params.dy);
                    
                    t_new[id] = fields.temperature[id] + params.dt * 
                        params.viscosity * laplacian_t;
                    
                    // Simple pressure correction (Poisson-like)
                    double div_u = 
                        (fields.velocity_u[fields.idx(i+1, j, k)] - fields.velocity_u[fields.idx(i-1, j, k)])/(2*params.dx) +
                        (fields.velocity_v[fields.idx(i, j+1, k)] - fields.velocity_v[fields.idx(i, j-1, k)])/(2*params.dy);
                    
                    p_new[id] = fields.pressure[id] - 0.1 * div_u;
                }
            }
        }
        
        // Apply boundary conditions (simple reflection)
        apply_boundary_conditions(u_new, v_new, p_new, t_new);
        
        // Update fields
        fields.velocity_u = u_new;
        fields.velocity_v = v_new;
        fields.pressure = p_new;
        fields.temperature = t_new;
        fields.compute_vorticity();
        
        current_iteration++;
        current_time += params.dt;
    }
    
    void apply_boundary_conditions(std::vector<double>& u, std::vector<double>& v,
                                   std::vector<double>& p, std::vector<double>& t) {
        // Simple no-slip boundary conditions on walls
        for (int k = 0; k < params.nz; k++) {
            for (int j = 0; j < params.ny; j++) {
                // Left wall
                int left = fields.idx(0, j, k);
                u[left] = 0.0;
                v[left] = 0.0;
                
                // Right wall
                int right = fields.idx(params.nx-1, j, k);
                u[right] = 0.0;
                v[right] = 0.0;
            }
        }
        
        for (int k = 0; k < params.nz; k++) {
            for (int i = 0; i < params.nx; i++) {
                // Bottom wall
                int bottom = fields.idx(i, 0, k);
                u[bottom] = 0.0;
                v[bottom] = 0.0;
                
                // Top wall (moving lid)
                int top = fields.idx(i, params.ny-1, k);
                u[top] = 1.0;  // Moving lid
                v[top] = 0.0;
            }
        }
    }
    
    Node create_conduit_mesh(int iteration) const {
        Node mesh;
        int n_total = params.nx * params.ny * params.nz;
        
        // ========== COORDINATE SET ==========
        mesh["coordsets/coords/type"] = "uniform";
        mesh["coordsets/coords/dims/i"] = params.nx;
        mesh["coordsets/coords/dims/j"] = params.ny;
        mesh["coordsets/coords/dims/k"] = params.nz;
        
        mesh["coordsets/coords/origin/x"] = params.x_min;
        mesh["coordsets/coords/origin/y"] = params.y_min;
        mesh["coordsets/coords/origin/z"] = params.z_min;
        
        mesh["coordsets/coords/spacing/dx"] = params.dx;
        mesh["coordsets/coords/spacing/dy"] = params.dy;
        mesh["coordsets/coords/spacing/dz"] = params.dz;
        
        // ========== TOPOLOGY ==========
        mesh["topologies/mesh/type"] = "uniform";
        mesh["topologies/mesh/coordset"] = "coords";
        
        // ========== FIELDS ==========
        // Velocity vector field
        mesh["fields/velocity/association"] = "vertex";
        mesh["fields/velocity/topology"] = "mesh";
        mesh["fields/velocity/values/u"].set(fields.velocity_u.data(), n_total);
        mesh["fields/velocity/values/v"].set(fields.velocity_v.data(), n_total);
        mesh["fields/velocity/values/w"].set(fields.velocity_w.data(), n_total);
        mesh["fields/velocity/units"] = "m/s";
        
        // Pressure scalar field
        mesh["fields/pressure/association"] = "vertex";
        mesh["fields/pressure/topology"] = "mesh";
        mesh["fields/pressure/values"].set(fields.pressure.data(), n_total);
        mesh["fields/pressure/units"] = "Pa";
        
        // Temperature scalar field
        mesh["fields/temperature/association"] = "vertex";
        mesh["fields/temperature/topology"] = "mesh";
        mesh["fields/temperature/values"].set(fields.temperature.data(), n_total);
        mesh["fields/temperature/units"] = "K";
        
        // Vorticity scalar field
        mesh["fields/vorticity/association"] = "vertex";
        mesh["fields/vorticity/topology"] = "mesh";
        mesh["fields/vorticity/values"].set(fields.vorticity.data(), n_total);
        mesh["fields/vorticity/units"] = "1/s";
        
        // Velocity magnitude (derived field)
        mesh["fields/velocity_magnitude/association"] = "vertex";
        mesh["fields/velocity_magnitude/topology"] = "mesh";
        mesh["fields/velocity_magnitude/values"].set(DataType::float64(n_total));
        float64_array vel_mag = mesh["fields/velocity_magnitude/values"].value();
        
        for (int i = 0; i < n_total; i++) {
            vel_mag[i] = sqrt(
                fields.velocity_u[i] * fields.velocity_u[i] +
                fields.velocity_v[i] * fields.velocity_v[i] +
                fields.velocity_w[i] * fields.velocity_w[i]
            );
        }
        mesh["fields/velocity_magnitude/units"] = "m/s";
        
        // ========== STATE INFORMATION ==========
        mesh["state/cycle"] = iteration;
        mesh["state/time"] = current_time;
        mesh["state/iteration"] = iteration;
        mesh["state/dt"] = params.dt;
        
        // ========== SIMULATION METADATA ==========
        mesh["info/nx"] = params.nx;
        mesh["info/ny"] = params.ny;
        mesh["info/nz"] = params.nz;
        mesh["info/dx"] = params.dx;
        mesh["info/dy"] = params.dy;
        mesh["info/dz"] = params.dz;
        mesh["info/viscosity"] = params.viscosity;
        
        // ========== FIELD STATISTICS ==========
        compute_field_statistics(mesh);
        
        return mesh;
    }
    
    void compute_field_statistics(Node& mesh) const {
        int n_total = params.nx * params.ny * params.nz;
        
        // Pressure statistics
        double p_min = fields.pressure[0];
        double p_max = fields.pressure[0];
        double p_sum = 0.0;
        
        for (double p : fields.pressure) {
            p_min = std::min(p_min, p);
            p_max = std::max(p_max, p);
            p_sum += p;
        }
        
        mesh["statistics/pressure/min"] = p_min;
        mesh["statistics/pressure/max"] = p_max;
        mesh["statistics/pressure/mean"] = p_sum / n_total;
        
        // Velocity magnitude statistics
        double v_min = 1e100;
        double v_max = -1e100;
        double v_sum = 0.0;
        
        for (int i = 0; i < n_total; i++) {
            double v = sqrt(
                fields.velocity_u[i] * fields.velocity_u[i] +
                fields.velocity_v[i] * fields.velocity_v[i]
            );
            v_min = std::min(v_min, v);
            v_max = std::max(v_max, v);
            v_sum += v;
        }
        
        mesh["statistics/velocity_magnitude/min"] = v_min;
        mesh["statistics/velocity_magnitude/max"] = v_max;
        mesh["statistics/velocity_magnitude/mean"] = v_sum / n_total;
        
        // Kinetic energy
        double kinetic_energy = 0.0;
        for (int i = 0; i < n_total; i++) {
            kinetic_energy += 0.5 * (
                fields.velocity_u[i] * fields.velocity_u[i] +
                fields.velocity_v[i] * fields.velocity_v[i]
            );
        }
        mesh["statistics/kinetic_energy"] = kinetic_energy;
    }
    
    int get_current_iteration() const { return current_iteration; }
    double get_current_time() const { return current_time; }
    const SimulationFields& get_fields() const { return fields; }
};

// ==============================================
// DATA SAVER (USING CONDUIT)
// ==============================================
class SimulationDataSaver {
private:
    std::string output_dir;
    Node simulation_root;  // Root node containing all data
    
public:
    SimulationDataSaver(const std::string& dir) : output_dir(dir) {
        // Create output directory
        system(("mkdir -p " + output_dir).c_str());
        
        // Initialize simulation metadata
        simulation_root["metadata/simulation_type"] = "CFD_Demo";
        simulation_root["metadata/date"] = "2024-01-20";
        simulation_root["metadata/author"] = "CFD_Simulator";
        simulation_root["metadata/description"] = "10-iteration structured mesh simulation";
    }
    
    void save_iteration(const Node& mesh, int iteration, double time) {
        std::string iteration_key = "iterations/iteration_" + 
                                   std::to_string(iteration);
        
        // Store complete mesh data
        simulation_root[iteration_key + "/mesh"] = mesh;
        simulation_root[iteration_key + "/time"] = time;
        
        // Also save individual file for this iteration
        std::string filename = output_dir + "/iteration_" + 
                              std::to_string(iteration) + ".yaml";
        relay::io::save(mesh, filename, "yaml");
        
        // Save to HDF5 as well (binary, more efficient)
        // std::string h5_filename = output_dir + "/iteration_" + 
        //                          std::to_string(iteration) + ".hdf5";
        // relay::io::save(mesh, h5_filename);
        
        std::cout << "  Saved iteration " << iteration 
                  << " to " << filename << std::endl;
    }
    
    void save_summary() {
        // Save complete simulation data
        std::string summary_file = output_dir + "/simulation_summary.yaml";
        relay::io::save(simulation_root, summary_file, "yaml");
        
        // Save to HDF5
        // std::string h5_file = output_dir + "/simulation_summary.hdf5";
        // relay::io::save(simulation_root, h5_file);
        
        std::cout << "\nComplete simulation data saved to:" << std::endl;
        std::cout << "  YAML: " << summary_file << std::endl;
        // std::cout << "  HDF5: " << h5_file << std::endl;
    }
    
    void save_checkpoint(const SimulationEngine& engine, 
                         const std::string& checkpoint_name) {
        Node checkpoint;
        
        // Save simulation state
        checkpoint["state/iteration"] = engine.get_current_iteration();
        checkpoint["state/time"] = engine.get_current_time();
        
        // Save parameters
        checkpoint["params/nx"] = 20;
        checkpoint["params/ny"] = 15;
        checkpoint["params/nz"] = 10;
        checkpoint["params/dx"] = 0.1;
        checkpoint["params/dt"] = 0.01;
        
        // Save fields
        const SimulationFields& fields = engine.get_fields();
        int n_total = 20 * 15 * 10;  // Hardcoded for example
        
        checkpoint["fields/velocity_u"].set(fields.velocity_u.data(), n_total);
        checkpoint["fields/velocity_v"].set(fields.velocity_v.data(), n_total);
        checkpoint["fields/pressure"].set(fields.pressure.data(), n_total);
        checkpoint["fields/temperature"].set(fields.temperature.data(), n_total);
        
        // Save checkpoint
        // std::string checkpoint_file = output_dir + "/" + checkpoint_name + ".hdf5";
        // relay::io::save(checkpoint, checkpoint_file);
        
        // std::cout << "Checkpoint saved: " << checkpoint_file << std::endl;
    }
};

// ==============================================
// POST-PROCESSING UTILITIES
// ==============================================
class PostProcessor {
public:
    static void extract_slice(const Node& mesh, int k_slice, 
                              const std::string& output_file) {
        Node slice;
        
        // Extract 2D slice from 3D mesh
        int nx = mesh["info/nx"].to_int();
        int ny = mesh["info/ny"].to_int();
        int n_slice = nx * ny;
        
        // Create 2D coordinate set
        slice["coordsets/coords/type"] = "uniform";
        slice["coordsets/coords/dims/i"] = nx;
        slice["coordsets/coords/dims/j"] = ny;
        slice["coordsets/coords/origin/x"] = mesh["coordsets/coords/origin/x"].to_double();
        slice["coordsets/coords/origin/y"] = mesh["coordsets/coords/origin/y"].to_double();
        slice["coordsets/coords/spacing/dx"] = mesh["coordsets/coords/spacing/dx"].to_double();
        slice["coordsets/coords/spacing/dy"] = mesh["coordsets/coords/spacing/dy"].to_double();
        
        slice["topologies/slice/type"] = "uniform";
        slice["topologies/slice/coordset"] = "coords";
        
        // Extract fields at this slice
        const float64_array& pressure = mesh["fields/pressure/values"].value();
        const float64_array& vel_u = mesh["fields/velocity/values/u"].value();
        const float64_array& vel_v = mesh["fields/velocity/values/v"].value();
        
        slice["fields/pressure/association"] = "vertex";
        slice["fields/pressure/topology"] = "slice";
        slice["fields/pressure/values"].set(DataType::float64(n_slice));
        
        slice["fields/velocity_u/association"] = "vertex";
        slice["fields/velocity_u/topology"] = "slice";
        slice["fields/velocity_u/values"].set(DataType::float64(n_slice));
        
        slice["fields/velocity_v/association"] = "vertex";
        slice["fields/velocity_v/topology"] = "slice";
        slice["fields/velocity_v/values"].set(DataType::float64(n_slice));
        
        float64_array p_slice = slice["fields/pressure/values"].value();
        float64_array u_slice = slice["fields/velocity_u/values"].value();
        float64_array v_slice = slice["fields/velocity_v/values"].value();
        
        // Copy data for this z-slice
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int src_idx = i + j * nx + k_slice * nx * ny;
                int dst_idx = i + j * nx;
                
                p_slice[dst_idx] = pressure[src_idx];
                u_slice[dst_idx] = vel_u[src_idx];
                v_slice[dst_idx] = vel_v[src_idx];
            }
        }
        
        relay::io::save(slice, output_file, "yaml");
        std::cout << "  Saved slice at k=" << k_slice 
                  << " to " << output_file << std::endl;
    }
    
    static void compute_statistics_over_time(const Node& simulation_root) {
        std::cout << "\n=== TIME SERIES STATISTICS ===" << std::endl;
        
        NodeConstIterator itr = simulation_root["iterations"].children();
        while (itr.has_next()) {
            const Node& iteration = itr.next();
            int cycle = iteration["mesh/state/cycle"].to_int();
            double time = iteration["mesh/state/time"].to_double();
            
            double kinetic_energy = iteration["mesh/statistics/kinetic_energy"].to_double();
            double max_vel = iteration["mesh/statistics/velocity_magnitude/max"].to_double();
            
            std::cout << std::fixed << std::setprecision(4)
                      << "Iter " << std::setw(3) << cycle 
                      << " | Time: " << std::setw(6) << time
                      << " | KE: " << std::setw(10) << kinetic_energy
                      << " | Max Vel: " << std::setw(8) << max_vel
                      << std::endl;
        }
    }
};

// ==============================================
// MAIN SIMULATION LOOP
// ==============================================
int main() {
    std::cout << "==========================================" << std::endl;
    std::cout << "STRUCTURED MESH SIMULATION - 10 ITERATIONS" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    // Set simulation parameters
    SimulationParams params;
    params.nx = 20;
    params.ny = 15;
    params.nz = 10;
    params.total_iterations = 10;
    
    // Create simulation components
    SimulationEngine engine(params);
    SimulationDataSaver saver("simulation_output");
    
    // Run simulation for 10 iterations
    for (int iter = 0; iter < params.total_iterations; iter++) {
        std::cout << "\n--- Iteration " << (iter + 1) << " ---" << std::endl;
        
        // Run one iteration
        engine.run_iteration();
        
        // Get current state
        int current_iter = engine.get_current_iteration();
        double current_time = engine.get_current_time();
        
        // Create Conduit mesh representation
        Node mesh = engine.create_conduit_mesh(current_iter);
        
        // Save iteration data
        saver.save_iteration(mesh, current_iter, current_time);
        
        // Save a slice for visualization (every 2 iterations)
        if (iter % 2 == 0) {
            std::string slice_file = "simulation_output/slice_iter_" + 
                                    std::to_string(current_iter) + "_k5.yaml";
            PostProcessor::extract_slice(mesh, 5, slice_file);
        }
        
        // Display some statistics
        std::cout << "  Time: " << current_time 
                  << ", Cycle: " << current_iter << std::endl;
    }
    
    // Save checkpoint
    saver.save_checkpoint(engine, "final_checkpoint");
    
    // Save complete simulation summary
    saver.save_summary();
    
    // Load and analyze saved data
    std::cout << "\n=== LOADING AND VERIFYING SAVED DATA ===" << std::endl;
    
    Node loaded_data;
    relay::io::load("simulation_output/iteration_5.yaml", "yaml", loaded_data);
    
    if (blueprint::mesh::verify(loaded_data, loaded_data)) {
        std::cout << "Loaded mesh verified successfully!" << std::endl;
        
        // Display some info
        std::cout << "Mesh dimensions: " 
                  << loaded_data["coordsets/coords/dims/i"].to_int() << " x "
                  << loaded_data["coordsets/coords/dims/j"].to_int() << " x "
                  << loaded_data["coordsets/coords/dims/k"].to_int() << std::endl;
        
        std::cout << "Number of fields: " 
                  << loaded_data["fields"].number_of_children() << std::endl;
    }
    
    // Load complete simulation summary
    // Node simulation_summary;
    // relay::io::load("simulation_output/simulation_summary.hdf5", simulation_summary);
    
    // Compute time series statistics
    // PostProcessor::compute_statistics_over_time(simulation_summary);
    
    std::cout << "\n==========================================" << std::endl;
    std::cout << "SIMULATION COMPLETED SUCCESSFULLY!" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    return 0;
}