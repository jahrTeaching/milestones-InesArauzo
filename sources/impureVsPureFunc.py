#Una función pura es una funcion que 
# SE PUUEDE SIMPLIFICAR g(x) + g(x)*(f(x)-1) =g(x)*h(x)
# NO TIENE EFECTOS COLATERALES, NO INTERACTUA CON EL MUNDO EXTERIOR (NO GLOBAL)

counter = 0

def impure1 (x):
    global counter
    counter = counter + 1 # esto es = que: counter += 1

    return x + counter

print('impure1(10)=  ', impure1(10))
print('impure1(10)=  ', impure1(10))

#PIDO 2 VECES LO MISMO Y EL RESULTADO ES DISTINTO PQ LA FUNCION ES IMPURA Y ME HA MODIFICADO EL GLOBAL

def impure2(l):
    l += 1
    print('id_dentro',id(l))
    return l

print('impure2(10)=  ', impure2(10))
print('impure2(10)=  ', impure2(10))


def pure(l):
    l = l + 2
    return l

print('pure(10) = ',  pure(10))


def Fibonacci_recurssive(n): # esto no es efectivo prque va haciendo una pila de llamadas y almacenando los resultados NO RENTA, PERO SE PUEDE PALIAR
    return Fibonacci(n - 1) + Fibonacci(n - 2)
