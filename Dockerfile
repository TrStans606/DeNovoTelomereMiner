# Use the conda base image
FROM conda/miniconda3

# Set the working directory
WORKDIR /docker

# Copy the environment file
COPY environment.yml .

# Create the environment
RUN conda env create -f environment.yml

# Activate the environment
RUN echo "source activate de_novo" > ~/.bashrc
ENV PATH /opt/conda/envs/de_novo/bin:$PATH

# Copy the repository content
COPY . .

# Clone the repository
CMD ["git", "clone", "https://github.com/TrStans606/DeNovoTelomereMiner"]
