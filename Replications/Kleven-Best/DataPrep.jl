module DataPrep

export prims # Export primitives from the data

# Load libraries
using Parameters
using CSV
using DataFrames
using Plots
using Trapz
using NumericalIntegration
using SpecialFunctions

# Load data
dta = DataFrame(CSV.File("data/statae050.csv", header = false));

ndta = size(dta)[1]; # Data size

# Rename columns
rename!(dta, [:z, :Hzy, :hzy, :Hzo, :hzo, :mtrz, :ω, :jω, :Jω, :dtapars]);

# Add a struct of parameters from the data
@with_kw struct Dta_params
    e
    d ::Float64 = 0 # Benchmark scenario where δ=0
    zymin
    αzy
    zomin
    αzo
    γ ::Float64 = 10
    R ::Int64 = 4000
end;

dtapars = Dta_params(e = dta.dtapars[1],
                     zymin = dta.dtapars[2],
                     αzy = dta.dtapars[3],
                     zomin = dta.dtapars[4],
                     αzo = dta.dtapars[5]);

# Function to make a distribution monotonic
function make_monotone(dist)
    new = copy(dist)
    ndata = length(dist)
    for i in 2:ndata
        if new[i] <= new[i-1]
            new[i] = new[i-1]
            nextind = findfirst(new .> new[i-1])[2]
            next = new[nextind]
            new[i] = (next + new[i-1]*(nextind-i)) / (nextind-i+1)
        end
    end
    return new
end;

# Function that smoothes distributions
# (Same procedure as in Saez)
function smooth_dist(dist, niter)
    old = copy(dist)
    new = copy(dist)
    ndata = length(dist)
    for i in 1:niter
        for j in 2:ndata-1
            new[j] = 0.3*old[j-1] + 0.4*old[j] + 0.3*old[j+1]
        end
        old = copy(new)
    end
    return new
end;

# Find the distribution of ω(n, z_y/n)
# Do this by matching the CDF of ω and the CDF of z_y
# (Minimize the distance bewteen the CDF values)

# Get the difference between the Hzy and Jω
ωzdiff = abs.(repeat(dta.Hzy', ndta, 1) - repeat(dta.Jω, 1, ndta));

# Find the minimum difference indices for each Hzy (columnwise)
# ωinds[i] = min_j |Hzy[i] - J(ω)[j])|
_, ωinds = findmin(ωzdiff, dims = 1);
ωinds = getindex.(ωinds[1,:], 1)';

# Create ω distributions of these minimum differences
ωz = dta.ω[ωinds];
Jωzy = dta.Jω[ωinds];
jωzy = dta.jω[ωinds];

# Find the distribution of z_o(z_y)
# Follow a similar procedure as for finding the distribution of ω,
# this time matching the CDFs of z_y and z_o

# Difference between young z CDF and old z CDF
zozdiff = abs.(repeat(dta.Hzy', ndta, 1) - repeat(dta.Hzo, 1, ndta));

# Find the minimum difference indices for each Hzy (columnwise)
# zoinds[i] = min_j |Hzy[i] - Hzo[j])|
_, zoinds = findmin(zozdiff, dims = 1);
zoinds = getindex.(zoinds[1,:], 1)';

# Create z distributions of these minimum differences
zoz = dta.z[zoinds];
Hzoz = dta.Hzo[zoinds];
hzoz = dta.hzo[zoinds];

# Replace z_o(z_y) for top incomes z>150k

# Split the z vector at 150,000
znew1 = dta.z[dta.z.<150000];
znew2 = dta.z[dta.z.>150000];
nz1 = size(znew1)[1];

# Construct z_o(z_y) incomes above 150,000
zohigh = dtapars.zomin * (znew2 / dtapars.zymin).^(dtapars.αzy/dtapars.αzo);

# Append to incomes below 150,000
zoznew = [zoz[1:nz1]; zohigh];

# Compute ω_0, make ω_0 monotonic, and smooth ω_0
ω0 = copy(ωz);
ω0 = make_monotone(ωz);
ω0 = smooth_dist(ω0, 2000);

# Compute ω_0, make ω_0 monotonic, and smooth ω_0
ω0 = copy(ωz);
ω0 = make_monotone(ωz);
ω0 = smooth_dist(ω0, 2000);

# Compute the ability levels, make monotone, and smooth
# (n = ability)
τ = dta.mtrz./100;
n = (dta.z.^(1/(1+dtapars.e))) .* ((1 .- τ).^(-dtapars.e/(1+dtapars.e)));
n = make_monotone(n');
n = smooth_dist(n', 2000);


# Compute the ability distribution and smooth
fn = diff(dta.Hzy) ./ diff(n);
push!(fn, 0);
fn = smooth_dist(fn, 3000);

# Create the CDF and normalize so sum = 1
Fn = cumul_integrate(n, fn); # CDF
fn = fn/Fn[end]; # Normalize
Fn = cumul_integrate(n, fn); # Normalize

# Pareto tails for the ability distribution
# Above z = 2,000,000

# Calculate the pareto distribution
haz = (n.*fn)./(1 .- Fn); # Clarify what this is?
pareto_ind = findfirst(dta.z .> 2000000);
pareto_α = haz[pareto_ind];
φ = 1 - Fn[pareto_ind];
pareto_lb = n[pareto_ind] * φ^(1/pareto_α);

# Update the distributions fn, Fn
fn[pareto_ind:end] = pareto_α .* pareto_lb^pareto_α ./ n[pareto_ind:end].^(1+pareto_α);
Fn[pareto_ind:end] = 1 .- (pareto_lb^pareto_α ./ n[pareto_ind:end].^pareto_α);

# Trim the ability distribution at 27,200
ntop = findfirst(n .> 27200) - 1;
n = n[1:ntop];
fnnew = fn[1:ntop];
Fnnew = Fn[1:ntop];
ω0 = ω0[1:ntop];

# Re-normalize distributions
fn = fnnew / Fnnew[ntop];
Fn = Fnnew / Fnnew[ntop];

# Squeeze the top of the z grid (clarify?)
z = zeros(ntop, 1);
for i = 1:ntop
    z[i] = exp(5 + ((i-1) * (11.52/(ntop-1))))
end;

# Define η (from hamiltonian.agedep.m, easier to do here)
ω0prime = diff(ω0) ./ diff(n);
push!(ω0prime, ω0prime[end]);
η = dtapars.d .+ (ω0prime .* n ./ ω0);

# Create a struct of data primitives to return (when modularized)
@with_kw struct Prims_struct
    z
    n
    ω0
    fn
    ntop
    γ
    d
    e
    η
end;

prims = Prims_struct(z, n, ω0, fn, ntop, dtapars.γ, dtapars.d, dtapars.e, η)

end