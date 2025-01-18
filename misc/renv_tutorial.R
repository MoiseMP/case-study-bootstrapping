# Using Renv

#install.packages('renv')

# Initialize
renv::init()

# Check Status
renv::status()

# View active library paths
.libPaths()


# To exit renv
# renv::deactivate()

# To remove renv
renv::remove()