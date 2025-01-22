# Using Renv

#install.packages('renv')

# Initialize
renv::init()

# Check Status
renv::status()

# View active library paths
.libPaths()

# Snapshot current state
renv::snapshot()

# Make environment in sync with Renv
renv::restore()

# To exit renv
# renv::deactivate()

# To remove renv
renv::remove()