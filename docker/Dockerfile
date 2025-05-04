# Use the conda base image
FROM continuumio/miniconda3:main

# Set the working directory
WORKDIR /docker

# Install Git
RUN apt update && apt upgrade && apt install -y git build-essential curl libboost-all-dev --no-install-recommends

# Clone the repository
RUN git clone https://github.com/TrStans606/DeNovoTelomereMiner 
WORKDIR /docker/DeNovoTelomereMiner

# Create the environment
RUN conda env create -f environment.yml

# Activate the environment
RUN echo "source activate de_novo" > ~/.bashrc
ENV PATH=/opt/conda/envs/de_novo/bin:$PATH

# Copy the repository content
COPY . .

#compile teloport
RUN bash install.sh

# Clone the repository
CMD ["bash", "-c", "source ~/.bashrc", "&&", "python3 DeNovoTelomereMiner.py -h"]
