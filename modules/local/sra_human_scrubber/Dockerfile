FROM ncbi/sra-human-scrubber:2.2.1

ENV DEBIAN_FRONTEND=noninteractive

# Install seqtk, pigz, and basic tools useful in NF processes
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        seqtk \
        pigz \
        bash \
        coreutils \
        procps \
        ca-certificates \
        curl && \
    rm -rf /var/lib/apt/lists/*