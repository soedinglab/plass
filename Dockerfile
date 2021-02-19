FROM debian:stable-slim as builder
RUN apt-get update && apt-get install -y \
      build-essential cmake xxd git zlib1g-dev libbz2-dev \
      && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/source
ADD . .

WORKDIR /opt/source/build_sse2
RUN cmake -DHAVE_MPI=0 -DHAVE_TESTS=0 -DHAVE_SSE2=1 -DCMAKE_BUILD_TYPE=Release ..
RUN make -j $(nproc --all)

WORKDIR /opt/source/build_sse41
RUN cmake -DHAVE_MPI=0 -DHAVE_TESTS=0 -DHAVE_SSE4_1=1 -DCMAKE_BUILD_TYPE=Release ..
RUN make -j $(nproc --all)

WORKDIR /opt/source/build_avx2
RUN cmake -DHAVE_MPI=0 -DHAVE_TESTS=0 -DHAVE_AVX2=1 -DCMAKE_BUILD_TYPE=Release ..
RUN make -j $(nproc --all)

FROM debian:stable-slim
MAINTAINER Milot Mirdita <milot@mirdita.de>
RUN apt-get update && apt-get install -y \
      gawk bash grep libstdc++6 libgomp1 libatomic1 zlib1g libbz2-1.0 wget tar \
      && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/source/build_sse2/src/plass /usr/local/bin/plass_sse2
COPY --from=builder /opt/source/build_sse41/src/plass /usr/local/bin/plass_sse41
COPY --from=builder /opt/source/build_avx2/src/plass /usr/local/bin/plass_avx2
RUN echo '#!/bin/sh\n\
FLAGS="$(grep -m 1 "^flags" /proc/cpuinfo)"\n\
case "${FLAGS}" in\n\
  *avx2*) exec /usr/local/bin/plass_avx2 "$@" ;;\n\
  *sse4_1*) exec /usr/local/bin/plass_sse41 "$@" ;;\n\
  *) exec /usr/local/bin/plass_sse2 "$@" ;;\n\
esac' >> /usr/local/bin/plass
RUN chmod +x /usr/local/bin/plass

VOLUME ["/app"]
WORKDIR /app
ENTRYPOINT ["plass"]
