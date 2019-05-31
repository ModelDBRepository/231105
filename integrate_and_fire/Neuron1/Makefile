# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
#CC=icpc
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
#CFLAGS=-c -Wall -std=gnu++11


#CFLAGS=-c -O2 -std=c++0x
CFLAGS=-c -O2 -std=c++11


#CFLAGS=-c -std=gnu++11
#CFLAGS=-c -std=gnu++0x
#-std=c++0x
all: net

net: Random.o SystemConstants.o PoissonGen.o NeuronMatrix.o ODE_file.o NeuronReader.o Connections.o IntegrationMethod.o WB_Neuron.o
	$(CC) Random.o SystemConstants.o PoissonGen.o NeuronMatrix.o ODE_file.o NeuronReader.o Connections.o IntegrationMethod.o WB_Neuron.o -o net.ex
#Connections.o IntegrationMethod.o netMatrix.o netReader.o ODE_file.o PoissonGen.o Random.o SystemConstants.o WB_net.o
#	$(CC) Connections.o IntegrationMethod.o netMatrix.o netReader.o ODE_file.o PoissonGen.o Random.o SystemConstants.o #WB_net.o -o net.exe

IntegrationMethod.o: IntegrationMethod.cpp Connections.o ODE_file.o
	$(CC) $(CFLAGS) IntegrationMethod.cpp

Connections.o: Connections.cpp NeuronMatrix.o
	$(CC) $(CFLAGS) Connections.cpp

NeuronMatrix.o: NeuronMatrix.cpp PoissonGen.o
	$(CC) $(CFLAGS) NeuronMatrix.cpp

NeuronReader.o: NeuronReader.cpp SystemConstants.o
	$(CC) $(CFLAGS) NeuronReader.cpp

ODE_file.o: ODE_file.cpp PoissonGen.o
	$(CC) $(CFLAGS) ODE_file.cpp

PoissonGen.o: PoissonGen.cpp Random.o SystemConstants.o
	$(CC) $(CFLAGS) PoissonGen.cpp

Random.o: Random.cpp gamma.o #added gamma 18.7.17
	$(CC) $(CFLAGS) Random.cpp
 
gamma.o: gamma.cpp
	$(CC) $(CFLAGS) gamma.cpp
 
SystemConstants.o: SystemConstants.cpp
	$(CC) $(CFLAGS) SystemConstants.cpp

WB_Neuron.o: WB_Neuron.cpp IntegrationMethod.o NeuronReader.o
	$(CC) $(CFLAGS) WB_Neuron.cpp

clean:
	rm *o net.ex

