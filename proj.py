import matplotlib.pyplot as plt

# używany do zrobienia wykresów

JacobiB = []
JacobiC = []
GaussSeidelB = []
GaussSeidelC = []

JacobiE = []
GaussSeidelE = []
LUE = []

iterationsJacobi = []
iterationsGaussSeidel = []
sizes = [100, 500, 1000, 2000, 3000, 4000, 5000]

#B

file = open("JacobiErrB.txt","r")

JacobiB = file.read().splitlines()

file.close()

for i in range(len(JacobiB)):
    JacobiB[i] = float(JacobiB[i])
    iterationsJacobi.append(i+1)

file = open("GaussSeidelB.txt","r")

GaussSeidelB = file.read().splitlines()

file.close()

for i in range(len(GaussSeidelB)):
    GaussSeidelB[i] = float(GaussSeidelB[i])
    iterationsGaussSeidel.append(i+1)
    
plt.semilogy(iterationsJacobi,JacobiB, label="Jacobi")
plt.semilogy(iterationsGaussSeidel,GaussSeidelB,label="Gauss-Seidel")
plt.xlabel("Iteracja")
plt.ylabel("Błąd")
plt.legend()
plt.show()    
    
# C
    
iterationsJacobi = []
iterationsGaussSeidel = []
    
file = open("JacobiErrC.txt","r")

JacobiC = file.read().splitlines()

file.close()

for i in range(len(JacobiC)):
    JacobiC[i] = float(JacobiC[i])
    iterationsJacobi.append(i+1)    

file = open("GaussSeidelC.txt","r")

GaussSeidelC = file.read().splitlines()

file.close()

for i in range(len(GaussSeidelC)):
    GaussSeidelC[i] = float(GaussSeidelC[i])
    iterationsGaussSeidel.append(i+1)
    
plt.semilogy(iterationsJacobi,JacobiC, label="Jacobi")
plt.semilogy(iterationsGaussSeidel,GaussSeidelC,label="Gauss-Seidel")
plt.xlabel("Iteracja")
plt.ylabel("Błąd")
plt.legend()
plt.show()  

# E

file = open("TimeE.txt","r")

for i in range(21):
    line = file.readline()
    line = line.rstrip('\n')
    if(i % 3 == 0):
        JacobiE.append(float(line)/ 1000)
    elif(i%3 == 1):
        GaussSeidelE.append(float(line)/ 1000)
    else:
        LUE.append(float(line) / 1000)
file.close()

plt.plot(sizes,JacobiE, label="Jacobi")
plt.plot(sizes,GaussSeidelE,label="Gauss-Seidel")
plt.xlabel("Rozmiar Macierzy")
plt.ylabel("Czas w sekundach")
plt.legend()
plt.show() 

plt.plot(sizes,JacobiE, label="Jacobi")
plt.plot(sizes,GaussSeidelE,label="Gauss-Seidel")
plt.plot(sizes,LUE, label="LU")
plt.xlabel("Rozmiar Macierzy")
plt.ylabel("Czas w sekundach")
plt.legend()
plt.show() 

