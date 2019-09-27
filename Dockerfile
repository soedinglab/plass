FROM debian:stable-slim as builder
RUN apt-get update && apt-get install -y \
      build-essential cmake xxd git zlib1g-dev libbz2-dev \
      && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/source
ADD . .

WORKDIR /opt/source/build_sse
RUN cmake -DHAVE_MPI=0 -DHAVE_TESTS=0 -DHAVE_SSE4_1=1 -DCMAKE_BUILD_TYPE=Release ..
RUN make -j $(nproc --all)

WORKDIR /opt/source/build_avx
RUN cmake -DHAVE_MPI=0 -DHAVE_TESTS=0 -DHAVE_AVX2=1 -DCMAKE_BUILD_TYPE=Release ..
RUN make -j $(nproc --all)

FROM debian:stable-slim
MAINTAINER Milot Mirdita <milot@mirdita.de>
RUN apt-get update && apt-get install -y \
      gawk bash grep libstdc++6 libgomp1 zlib1g libbz2-1.0 \
      && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/source/build_sse/src/plass /usr/local/bin/plass_sse42
COPY --from=builder /opt/source/build_avx/src/plass /usr/local/bin/plass_avx2
RUN echo '#!/bin/bash\n\
if $(grep -q -E "^flags.+avx2" /proc/cpuinfo); then\n\
    exec /usr/local/bin/plass_avx2 "$@"\n\
else\n\
    exec /usr/local/bin/plass_sse42 "$@"\n\
fi'\
>> /usr/local/bin/plass
RUN chmod +x /usr/local/bin/plass

VOLUME ["/app"]
WORKDIR /app
ENTRYPOINT ["plass"]
