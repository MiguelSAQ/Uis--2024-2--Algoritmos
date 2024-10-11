import numpy as np
import random as rnd
from time import perf_counter
import matplotlib.pyplot as plt

Hasta = int(input("Hasta que numero: ")) # Número máximo de elementos que va a tener el último vector
Numero = int(input("Numero de elemtos: ")) # La cantidad de vectores que se van a generar
Desde = Hasta//Numero # En qué numero empieza (sin el 0) range y linscape
Paso = Hasta//Numero # De cuánto en cuánto aumenta los elementos del vector.

# Listas vacias donde se pondran los tiempos que se demore cada algoritmo de ordenamiento.
t_Insertion = []
t_Selection = []
t_Quick = []
t_Bubble = []


# - - - Inserion Sort - - - #

def insertion(A):
    i = 1
    while i < len(L):
        j = i
        while j > 0 and L[j - 1] > L[j]:
            L[j], L[j - 1] = L[j - 1], L[j]
            j -= 1
        i += 1
    return L


# - - - Selection Sort - - - #

def selection(A):
    for i in range(len(L) - 1):
        min = i
        for j in range(i + 1, len(L)):
            if (L[min] > L[j]):
                min = j
        if min != i:
            L[min], L[i] = L[i], L[min]
    return L 


# - - - Quick Sort - - - #

def quick(A):
    def partition(array, low, high):
        m = int((low + high) / 2)
        if array[m] < array[low]: 
            array[low], array[m] = array[m], array[low]
        if array[high] < array[low]:
            array[low], array[high] = array[high], array[low]
        if array[m] < array[high]:
            array[m], array[high] = array[high], array[m] 
        pivot = array[high]
        
        i = low - 1
        
        for j in range(low, high):
            if array[j] <= pivot:

                i = i + 1
                array[i], array[j] = array[j], array[i]
                
        array[i + 1], array[high] = array[high], array[i + 1]

        return i + 1

    def quickSort(array, low, high):
        if low < high:

            pi = partition(array, low, high)
            quickSort(array, low, pi - 1)
            quickSort(array, pi + 1, high)

    quickSort(L, 0, len(L) - 1)
    return L


# - - - Bubble Sort - - - #

def bubble(A):
    for i in range(len(L) - 1):
        for j in range(len(L) - i - 1):
            if L[j] > L[j + 1]:
                L[j], L[j + 1] = L[j + 1], L[j]
    return L


# - - - Ordenamiento - - - #

K=1
for n in range(Desde, Hasta+1, Paso) : 

    # - - -  Random Vector - - - #

    L = []

    for i in range(n):
        L.append(int(100*rnd.random()))

    L_original = np.copy(L) # Guarda una copia de L en L_original

    # - - - Inserion Sort - - - #
            
    L = np.copy(L_original) # Define L como la copia del original
    ti1 = perf_counter()
    insertion(L)
    tf1 = perf_counter()  
    T1 = tf1 - ti1

    # - - - Selection Sort - - - #

    L = np.copy(L_original) # Define L como la copia del original
    ti2 = perf_counter()
    selection(L)   
    tf2 = perf_counter()
    T2 = tf2 - ti2

    # - - - Quick Sort - - - #

    L = np.copy(L_original) # Define L como la copia del original
    ti3 = perf_counter()
    quick(L)
    tf3 = perf_counter()
    T3 = tf3 - ti3

    # - - - Bubble sort - - - #

    L = np.copy(L_original) # Define L como la copia del original   
    ti4 = perf_counter()
    bubble(L)
    tf4 = perf_counter()
    T4 = tf4 - ti4

    # - - - Tiempos de ejecucion - - - #

    t_Insertion.append(T1)
    t_Selection.append(T2)
    t_Quick.append(T3)
    t_Bubble.append(T4)

    # - - - Numero del Paso - - - #

    print(f"\rPaso {K}/{Numero}", end='')
    K += 1

# - - - Tiempos de ejecucion - - - #

print([f"{x:.2f}" for x in t_Insertion])
print([f"{x:.2f}" for x in t_Selection])
print([f"{x:.2f}" for x in t_Quick])
print([f"{x:.2f}" for x in t_Bubble])

# - - - Grafíca del tiempo vs. el numero de datos - - - #

Datos = np.linspace(Desde, Hasta, Numero)

plt.figure(figsize=(10, 6))

# Ajuste cuadrático para cada conjunto de datos
coef_Insertion = np.polyfit(Datos, t_Insertion, 2)
coef_Selection = np.polyfit(Datos, t_Selection, 2)
coef_Quick = np.polyfit(Datos, t_Quick, 2)
coef_Bubble = np.polyfit(Datos, t_Bubble, 2)

print(f"Coeficientes Insertion: {coef_Insertion[0]}x^2 + {coef_Insertion[1]}x + {coef_Insertion[2]}")
print(f"Coeficientes Selection: {coef_Selection[0]}x^2 + {coef_Selection[1]}x +  {coef_Selection[2]}")
print(f"Coeficientes Quick: {coef_Quick[0]}x^2 + {coef_Quick[1]}x +  {coef_Quick[2]}")
print(f"Coeficientes Bubble: {coef_Bubble[0]}x^2 + {coef_Bubble[1]}x +  {coef_Bubble[2]}")


# Crear valores ajustados para cada conjunto de datos
ajuste_Insertion = np.polyval(coef_Insertion, Datos)
ajuste_Selection = np.polyval(coef_Selection, Datos)
ajuste_Quick = np.polyval(coef_Quick, Datos)
ajuste_Bubble = np.polyval(coef_Bubble, Datos)

# Graficar los datos originales
plt.scatter(Datos, t_Insertion, color='green', label='Insertion Sort')
plt.scatter(Datos, t_Selection, color='red', label='Selection Sort')
plt.scatter(Datos, t_Quick, color='cyan', label='Quick Sort')
plt.scatter(Datos, t_Bubble, color='blue', label='Bubble Sort')

# Graficar los ajustes cuadráticos
plt.plot(Datos, ajuste_Insertion, color='green', linestyle='--')
plt.plot(Datos, ajuste_Selection, color='red', linestyle='--')
plt.plot(Datos, ajuste_Quick, color='cyan', linestyle='--')
plt.plot(Datos, ajuste_Bubble, color='blue', linestyle='--')

plt.xlabel('Numero de datos a ordenar')
plt.ylabel('Tiempo en ordenar (s)')
plt.title('Comparación de Tiempos de Ejecución para Algoritmos de Ordenamiento​', fontweight='bold')
plt.legend()
plt.grid()
plt.show()