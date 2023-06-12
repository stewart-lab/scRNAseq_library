

# Set the base path and the name of the environment
conda_envs_path="/isiseqruns/jfreeman_tmp_home/bin/miniconda3/envs"
name="TEST_ENVIRONMENT"

# Construct the full path
full_path="${conda_envs_path}/${name}"

# Replace PREFIX_PLACEHOLDER with the desired path in the .yml file
sed "s|PREFIX_PLACEHOLDER|$full_path|" environment_template.yml > environment.yml

# Replace NAME_PLACEHOLDER with the desired name in the generated environment.yml file
sed -i "s|NAME_PLACEHOLDER|$name|" environment.yml

conda env create -f environment.yml

conda activate $name

