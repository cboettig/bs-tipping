
f=2
cJA=  1e-3
cJF=0.5 
cFA=0.3 
v=1 
h=8 
s=0.5
Fo=100
d = 0.1
qE = 1 
mA = 0


Fhat <- function(A) d*Fo / (d + cFA *A)

Jhat <- function(A) f*A / ( cJA * A + cJF* v * Fhat(A) / (h + v + cJF * Fhat(A))   )

Amotion <- function(A){ 
  - qE * A + 
    s * Jhat(A) - 
    mA * A 
  }



qE = 1
curve(Amotion, 0, 1100)

qE = 10
curve(Amotion, 0, 100)


qE = 20
curve(Amotion, 0, 100)
