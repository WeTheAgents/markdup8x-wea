# Picard-shim image: ships markdup-wea + a thin bash shim that mimics the
# Picard MarkDuplicates CLI for nf-core/modules drop-in testing.
# See: docs/nfcore-invocation.md, scripts/picard-shim.sh.

FROM rust:1.85-slim AS builder
WORKDIR /src
COPY Cargo.toml Cargo.lock ./
COPY src ./src
RUN cargo build --release --bin markdup-wea

FROM debian:stable-slim
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        samtools bash gawk procps ca-certificates \
 && rm -rf /var/lib/apt/lists/*
COPY --from=builder /src/target/release/markdup-wea /usr/local/bin/markdup-wea
COPY scripts/picard-shim.sh /usr/local/bin/picard
RUN chmod +x /usr/local/bin/picard /usr/local/bin/markdup-wea

ENV WEA_BIN=/usr/local/bin/markdup-wea
