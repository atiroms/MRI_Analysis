Sys.info()
R.Version()

# Memory usage information
memory.size(T)   # Maximum amount of memory obtained from the OS
memory.size(F)   # The amount currently in use
memory.size(NA)  # Memory limit

# Set memory limit
memory.limit(1000000)

# Test memory limit
x <- data.frame(x = runif(1000000000*4))

# Clear vector and memory
x<-NULL
gc()
