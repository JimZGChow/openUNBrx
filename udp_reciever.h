#ifndef UDP_RECIEVER_H
#define UDP_RECIEVER_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netinet/in.h>

class UDP_Reciever
{
public:
    UDP_Reciever(int maxDataSize);

    int readData(uint8_t* data);
private:
    int sock;
    struct sockaddr_in servaddr;
    int client_sock , c;
    struct sockaddr_in client;
    int maxDataSize;
};

#endif // UDP_RECIEVER_H
