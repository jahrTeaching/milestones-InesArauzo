from numpy import array, reshape


a = 4
pa = a

print(id(a), id(pa)) # los dos me dan el id de 3, estoy llamando al misnmo espacio de memoria, por eso 

a = 3 # el tres es un objeto

print(id(a), id(pa)) # como  el tres es un objeto distinto al 4 y yo la etiqueta pa se la estaba dando a 4, 
                     # ahora cambio el objetode pa pero no el de a, a a se le ha asignado un nuevo objeto


#############################################################################################################

L = [1, 2, 3] 
pL = L # alias

print("pL = ", pL)
print("idpL = ", id(L))
L[0] = 0

print("pL = ", pL)
print("idpL = ", id(L))
L = L + [1, 2, 3]
print("pL = ", L)
print("idL = ", id(L))

L = [1, 2, 3] 
pL = L[:] # cloning

#############################################################################################################

b = array([1, 2, 3] )
pb = b

print(id(b), id(pb))

b =b+1

print("b = ", b, "pb = ", pb)
print(id(b), id(pb))

b[:] =b[:]+1

print("b = ", b, "pb = ", pb)
print(id(b), id(pb))



#############################################################################################################

c = array([1, 2, 3])
pc = c[:].copy()

c[:] = c[:] + 1
print("c = ", c, "pc = ", pc)

#############################################################################################################

d = array([1, 2, 3, 4])

pd =reshape(d, (2, 2))

print("d = ", d, "pd = ", pd)
print(id(d), id(pd))

d[:] = 0

print("d = ", d, "pd = ", pd)
print(id(d), id(pd))

############################################

U = array([1, 2, 3, 4])
U = reshape(U, (2, 2))
print("U =", U)
#U[:,[0,1]] = U[:,[1,0]] #esto es un mapping y esuna forma nais de hacer un mapping
# U[:,0],  U[:,1] = U[:,1],  U[:,0]
U[:,0],  U[:,1] = U[:,1].copy(),  U[:,0].copy()
print("U =", U) 