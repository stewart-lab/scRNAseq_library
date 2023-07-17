# Change to the output directory
cd output

# Get the absolute path to the current directory
output_dir=$(pwd)

# Go back to the previous directory
cd ..

# Run the Docker container with the correct volume mapping
docker run -v ${output_dir}:/scRNA-seq/output -it jfreeman88/scrna-seq:latest
