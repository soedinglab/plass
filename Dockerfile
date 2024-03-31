FROM --platform=$BUILDPLATFORM debian:stable-slim as builder
ARG TARGETARCH

RUN dpkg --add-architecture $TARGETARCH \
    && apt-get update \
    && apt-get install -y \
      build-essential cmake xxd git \
      zlib1g-dev libbz2-dev libatomic1 \
      crossbuild-essential-$TARGETARCH zlib1g-dev:$TARGETARCH libbz2-dev:$TARGETARCH \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/build
ADD . .

RUN if [ "$TARGETARCH" = "arm64" ]; then \
      mkdir -p build_$TARGETARCH/src; \
      cd /opt/build/build_$TARGETARCH; \
      CC=aarch64-linux-gnu-gcc CXX=aarch64-linux-gnu-g++ cmake -DHAVE_ARM8=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
      make -j $(nproc --all); \
      mv src/plass /opt/build/plass_arch; \
      mv src/penguin /opt/build/penguin_arch; \
      touch /opt/build/plass_sse2 /opt/build/plass_sse41 /opt/build/plass_avx2; \
      touch /opt/build/penguin_sse2 /opt/build/penguin_sse41 /opt/build/penguin_avx2; \
    else \
      mkdir -p build_sse2/src && mkdir -p build_sse41/src && mkdir -p build_avx2/src; \
      cd /opt/build/build_sse2; \
      cmake -DHAVE_SSE2=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
      make -j $(nproc --all); \
      mv src/plass /opt/build/plass_sse2; \
      mv src/penguin /opt/build/penguin_sse2; \
      cd /opt/build/build_sse41; \
      cmake -DHAVE_SSE4_1=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
      make -j $(nproc --all); \
      mv src/plass /opt/build/plass_sse41; \
      mv src/penguin /opt/build/penguin_sse41; \
      cd /opt/build/build_avx2; \
      cmake -DHAVE_AVX2=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
      make -j $(nproc --all); \
      mv src/plass /opt/build/plass_avx2; \
      mv src/penguin /opt/build/penguin_avx2; \
      touch /opt/build/plass_arch; \
      touch /opt/build/penguin_arch; \
    fi

FROM debian:stable-slim
ARG TARGETARCH

RUN apt-get update && apt-get install -y \
      gawk bash grep libstdc++6 libgomp1 libatomic1 zlib1g libbz2-1.0 wget tar \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/build/plass_arch /opt/build/plass_sse2 /opt/build/plass_sse41 /opt/build/plass_avx2 /usr/local/bin/
COPY --from=builder /opt/build/penguin_arch /opt/build/penguin_sse2 /opt/build/penguin_sse41 /opt/build/penguin_avx2 /usr/local/bin/

RUN if [ "$TARGETARCH" != "arm64" ]; then \
      echo '#!/bin/sh\n\
      FLAGS="$(grep -m 1 "^flags" /proc/cpuinfo)"\n\
      case "${FLAGS}" in\n\
          *avx2*) exec /usr/local/bin/plass_avx2 "$@" ;;\n\
          *sse4_1*) exec /usr/local/bin/plass_sse41 "$@" ;;\n\
          *) exec /usr/local/bin/plass_sse2 "$@" ;;\n\
      esac' >> /usr/local/bin/plass; \
      chmod +x /usr/local/bin/plass; \
      sed 's/plass/penguin/g' /usr/local/bin/plass > /usr/local/bin/penguin; \
      chmod +x /usr/local/bin/penguin; \
    else \
      ln -s /usr/local/bin/plass_arch /usr/local/bin/plass; \
    fi
