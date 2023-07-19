# Use the official miniconda base image with Ubuntu 20.04
FROM continuumio/miniconda3:latest

# Update conda and install necessary system dependencies
RUN apt-get update -y && \
    apt-get install -y samtools && \
    apt-get install -y git && \
    apt-get install -y procps && \
    git clone https://github.com/cschlaffner/PoGo.git

# Set the working directory
WORKDIR /opt

# Copy the requirements.txt file and install Python dependencies
COPY requirements.txt /opt/requirements.txt
RUN pip install -r /opt/requirements.txt

# Install AGAT from the Bioconda channel
# RUN conda install -c bioconda agat

# Copy your Proteogenome script into the container
COPY Proteogenome.py /opt/Proteogenome.py

# Set the entry point for the Docker container
# ENTRYPOINT ["python", "proteogenome.py"]
