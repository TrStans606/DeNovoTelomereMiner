# Use the conda base image
FROM continuumio/miniconda3

# Set the working directory
WORKDIR /docker

# Install Git
RUN apt-get update && apt-get install -y git build-essential libboost-all-dev

# Copy the environment file
COPY environment.yml .

# Create the environment
RUN conda env create -f environment.yml

# Activate the environment
RUN echo "source activate de_novo" > ~/.bashrc
ENV PATH=/opt/conda/envs/de_novo/bin:$PATH

# Clone the repository
RUN git clone https://github.com/TrStans606/DeNovoTelomereMiner

# Copy the repository content
COPY . .

# Clone the repository
CMD ["bash", "-c", "source ~/.bashrc && git --version && bash install.sh && python3 DeNovoTelomereMiner.py -h"]
