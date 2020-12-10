#include "udp_reciever.h"

UDP_Reciever::UDP_Reciever(int _maxDataSize) {
/*
    std::cout << "Init UDP server (" << _maxDataSize << ")" << std::endl;
    // Creating socket file descriptor

    if ( (sock = socket(AF_INET, SOCK_DGRAM, 0)) < 0 ) {
        perror("socket creation failed");
        exit(EXIT_FAILURE);
    }

    memset(&servaddr, 0, sizeof(servaddr));
    memset(&client, 0, sizeof(client));

    // Filling server information
    servaddr.sin_family    = AF_INET; // IPv4
    servaddr.sin_addr.s_addr = INADDR_ANY;
    servaddr.sin_port = htons(12123);

    // Bind the socket with the server address
    if ( bind(sock, (const struct sockaddr *)&servaddr, sizeof(servaddr)) < 0 ) {
        perror("bind failed");
        exit(EXIT_FAILURE);
    }

    maxDataSize = _maxDataSize;
*/

/*
    std::cout << "Init TCP server (" << _maxDataSize << ")" << std::endl;
    //Create socket
    sock = socket(AF_INET , SOCK_STREAM , 0);
    if (sock == -1) {
        std::cout << "Could not create socket" << std::endl;
        return;
    }

    //Prepare the sockaddr_in structure
    servaddr.sin_family = AF_INET;
    servaddr.sin_addr.s_addr = INADDR_ANY;
    servaddr.sin_port = htons(12123);

    //Bind
    if( bind(sock,(struct sockaddr*)&servaddr , sizeof(servaddr)) < 0) {
        //print the error message
        std::cout << "bind failed. Error" << std::endl;
        perror("");
        return;
    }

    listen(sock , 3);

    //Accept and incoming connection
    std::cout << "Waiting for incoming connections..." << std::endl;
    c = sizeof(struct sockaddr_in);

    //accept connection from an incoming client
    client_sock = accept(sock, (struct sockaddr *)&client, (socklen_t*)&c);
    if (client_sock < 0) {
        std::cout << "accept failed" << std::endl;
        return;
    }
    std::cout << "Connection accepted" << std::endl;
*/
    std::cout << "Init TCP client (" << _maxDataSize << ")" << std::endl;
    sock = socket(AF_INET , SOCK_STREAM , 0);
    if (sock == -1) {
        std::cout << "Could not create socket" << std::endl;
        return;
    }

    //servaddr.sin_addr.s_addr = inet_addr("169.254.68.184");
    servaddr.sin_addr.s_addr = inet_addr("127.0.0.1");
    servaddr.sin_family = AF_INET;
    //servaddr.sin_addr.s_addr = INADDR_ANY;
    servaddr.sin_port = htons(12128);

    c = sizeof(servaddr);

    //accept connection from an incoming client
    client_sock = connect(sock, (sockaddr*)&servaddr, sizeof (servaddr));
    if (client_sock < 0) {
        std::cout << "accept failed" << std::endl;
        return;
    }
    std::cout << "Connection accepted" << std::endl;
}

int UDP_Reciever::readData(uint8_t *data) {


    int n = recv(sock, data, 1024*1024, MSG_WAITALL);
    //int n = recvfrom(sock, (char *)data, maxDataSize, MSG_WAITALL, ( struct sockaddr *) &client, (socklen_t*)&c);
    //std::cout << "Recieve " << n << std::endl;

    if (n < 0) {
        perror("err ");
        exit(-1);
    }

    n = n / (sizeof (float) * 2);

    //std::cout << n << std::endl;
    /*if (n == 4095) {
        std::cout << n << std::endl;
        ((float*)data)[4094 * 2] = 10;
        ((float*)data)[4094 * 2 + 1] = 10;
    }*/

    return n;
}
