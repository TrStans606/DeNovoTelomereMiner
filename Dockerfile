FROM conda/miniconda3

WORKDIR /docker

COPY environment.yml .

RUN conda env create -f environment.yml

RUN echo "source activate de_novo" > ~/.bashrc
ENV PATH /opt/conda/envs/de_novo/bin:$PATH

COPY . .

CMD ["git clone https://github.com/TrStans606/DeNovoTelomereMiner"]

