# Benchmarks

This is a transparent dense-design reference implementation, not an optimized
production estimator. Test fixture: 12 units x 10 periods; two sampled assignments
for cross-language inference. One BLAS/OpenMP thread; at most two Python processes.
No large Monte Carlo was launched. Formal 100k/500k performance benchmarks,
399/999-draw timing and memory profiling remain pending. No speed claims.
