FROM bioconductor/bioconductor_docker:RELEASE_3_20

ENV DEBIAN_FRONTEND=noninteractive
ENV HOST=0.0.0.0
ENV PORT=8788
ENV PYTHONUNBUFFERED=1

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        python3-venv \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /app/requirements.txt
RUN python3 -m pip install --break-system-packages --no-cache-dir -r /app/requirements.txt

RUN Rscript -e 'install.packages(c("dplyr","jsonlite","ggplot2","ggrepel","pheatmap","RColorBrewer","base64enc"), repos="https://cloud.r-project.org")' \
    && Rscript -e 'BiocManager::install(c("DESeq2","apeglm","limma","clusterProfiler","org.Hs.eg.db","org.Mm.eg.db","AnnotationDbi","enrichplot"), ask=FALSE, update=FALSE)'

COPY . /app

EXPOSE 8788

CMD ["python3", "server.py", "--host", "0.0.0.0"]
