FROM gcc:latest

WORKDIR /app

COPY matrix.h matrix.c test_matrix.c Makefile ./

RUN make test_matrix

CMD ["make", "test"]
