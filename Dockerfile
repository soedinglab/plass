FROM alpine:3.8 AS plass-builder
RUN apk add --no-cache gcc g++ cmake musl-dev vim git ninja zlib-dev bzip2-dev

WORKDIR /opt/plass
ADD . .

WORKDIR build_sse
RUN cmake -G Ninja -DHAVE_SSE4_1=1 -DCMAKE_BUILD_TYPE=Release ..
RUN ninja && ninja install

WORKDIR ../build_avx
RUN cmake -G Ninja -DHAVE_AVX2=1 -DCMAKE_BUILD_TYPE=Release ..
RUN ninja && ninja install

FROM alpine:3.8
MAINTAINER Milot Mirdita <milot@mirdita.de>
RUN apk add --no-cache gawk bash grep libstdc++ libgomp zlib libbz2

COPY --from=plass-builder /opt/plass/build_sse/src/plass /usr/local/bin/plass_sse42
COPY --from=plass-builder /opt/plass/build_avx/src/plass /usr/local/bin/plass_avx2
RUN echo -e '#!/bin/bash\n\
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
