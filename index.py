def mostrar_matriz(A):
    for fila in A:
        print(fila)

def ingresar_sistema():
    n = int(input("Numero de ecuaciones: "))

    if n <= 0:
        raise ValueError("El numero debe ser positivo")

    A = []
    b = []

    for i in range(n):
        fila = list(map(float, input(f"Coeficientes de la ecuacion {i+1}: ").split()))

        if len(fila) != n:
            raise ValueError("Cantidad incorrecta de coeficientes")

        A.append(fila)

    b = list(map(float, input("Terminos independientes: ").split()))

    if len(b) != n:
        raise ValueError("Cantidad incorrecta de terminos independientes")

    return A, b

def eliminacion_gaussiana(A, b):
    n = len(A)
    M = [A[i][:] + [b[i]] for i in range(n)]

    for i in range(n):
        p = max(range(i, n), key=lambda k: abs(M[k][i]))

        if abs(M[p][i]) < 1e-12:
            raise ValueError("Sistema sin solucion unica")

        M[i], M[p] = M[p], M[i]

        for j in range(i + 1, n):
            factor = M[j][i] / M[i][i]

            for k in range(i, n + 1):
                M[j][k] -= factor * M[i][k]

    x = [0.0] * n

    for i in range(n - 1, -1, -1):
        x[i] = (M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))) / M[i][i]

    return x

def gauss_jordan(A, b):
    n = len(A)
    M = [A[i][:] + [b[i]] for i in range(n)]

    for i in range(n):
        p = max(range(i, n), key=lambda k: abs(M[k][i]))

        if abs(M[p][i]) < 1e-12:
            raise ValueError("Sistema sin solucion unica")

        M[i], M[p] = M[p], M[i]

        pivote = M[i][i]

        for j in range(n + 1):
            M[i][j] /= pivote

        for k in range(n):
            if k != i:
                factor = M[k][i]

                for j in range(n + 1):
                    M[k][j] -= factor * M[i][j]

    return [M[i][n] for i in range(n)]

def gauss_seidel(A, b, tolerancia=1e-6, max_iter=100):
    n = len(A)

    if any(abs(A[i][i]) < 1e-12 for i in range(n)):
        raise ValueError("La diagonal contiene ceros")

    x = [0.0] * n

    for _ in range(max_iter):
        anterior = x[:]

        for i in range(n):
            suma = sum(A[i][j] * x[j] for j in range(n) if j != i)
            x[i] = (b[i] - suma) / A[i][i]

        error = max(abs(x[i] - anterior[i]) for i in range(n))

        if error < tolerancia:
            return x

    raise ValueError("El metodo no converge")


        if error < tolerancia:
            return x

    raise ValueError("El metodo no converge")

def ejecutar_pruebas():
    pruebas_directos = [
        ([[2, 1], [1, 3]], [5, 6]),
        ([[3, 2], [1, 2]], [16, 10]),
        ([[2, 1, 1], [0, 3, 2], [0, 0, 4]], [7, 12, 12])
    ]

    pruebas_iterativos = [
        ([[4, 1], [2, 5]], [9, 12]),
        ([[10, 1], [2, 8]], [21, 12]),
        ([[10, 1, 1], [1, 10, 1], [1, 1, 10]], [15, 24, 33])
    ]

    pruebas_incorrectas = [
        ([[0, 2], [1, 3]], [4, 5]),
        ([[1, 2], [2, 1]], [3, 3])
    ]

    metodos = [
        ("Eliminacion Gaussiana", eliminacion_gaussiana),
        ("Gauss-Jordan", gauss_jordan),
        ("Gauss-Seidel", gauss_seidel),
        ("Jacobi", jacobi)
    ]

    for nombre, metodo in metodos:
        print("\n", nombre)

        pruebas = pruebas_directos if metodo in [eliminacion_gaussiana, gauss_jordan] else pruebas_iterativos

        for i, (A, b) in enumerate(pruebas, 1):
            try:
                print(f"Correcto {i}:", metodo(A, b))
            except ValueError as e:
                print(f"Error:", e)

        for i, (A, b) in enumerate(pruebas_incorrectas, 1):
            try:
                print(f"Incorrecto {i}:", metodo(A, b))
            except ValueError as e:
                print(f"Error detectado {i}:", e)

def main():
    while True:
        print("\n===== METODOS NUMERICOS =====")
        print("1. Eliminacion Gaussiana")
        print("2. Gauss-Jordan")
        print("3. Gauss-Seidel")
        print("5. Ejecutar pruebas")
        print("6. Salir")

        opcion = input("Selecciona una opcion: ")

        if opcion == "6":
            print("Programa finalizado")
            break

        if opcion == "5":
            ejecutar_pruebas()
            continue

        metodos = {
            "1": eliminacion_gaussiana,
            "2": gauss_jordan,
            "3": gauss_seidel,
        }

        if opcion not in metodos:
            print("Opcion incorrecta")
            continue

        try:
            A, b = ingresar_sistema()
            resultado = metodos[opcion](A, b)

            print("\nSolucion del sistema:")
            for i, valor in enumerate(resultado):
                print(f"x{i+1} = {valor:.6f}")

        except (ValueError, IndexError, TypeError) as e:
            print("Error:", e)

if __name__ == "__main__":
    main()
