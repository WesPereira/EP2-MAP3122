################################################################################
##                          Métodos Numéricos e Aplicações                    ##
##                              Exercício Computacional                       ##
##                                                                            ##
##   Autores:                                                                 ##
##   Aldomar Pietro Silva Santana  NUSP: 10770162                             ##
##   Wesley Pereira de Almeida     NUSP: 10770093                             ##
##   NUSP : 10770093                                                          ##
##   Data : 16/02/2020                                                        ##
################################################################################
try:
    import numpy
    import math
    import matplotlib.pyplot as plt
except:
    print("[ERROR] Biblioteca não encontrada.")

def main():
    '''(None) -> (None)

    Função principal que comunica e mostra os resultados
    para o usuário.
    '''

    print("\t\tMétodo de Newton para calcular raízes de polinômios")
    print("\t\t\t  Exercício Computacional")
    print("\t\t\tMAP3122 - Quadimestral 2020")
    print("\t\t\t   Prof. Antoine Laurain")
    print("\t\t   Wesley P de Almeida  - NUSP: 10770093")
    print("\t\t   Aldomar P S Santana  - NUSP: 10770162")
    print()

    while True:
        print("Caso queira usar os valores padrões basta apertar ENTER.")

        #Definição das constantes
        c2 = 20
        nx = 100
        delta_x = 1/nx
        x_c = 0.7
        nt = 1000
        T = 1
        beta = 10
        delta_t = T/nt
        alpha = math.sqrt(c2)*(delta_t/delta_x)

        ti = 0.5
        tf = 1
        xk = [0.2, 0.3, 0.9]
        xr = 0.7
        K = 3
        #Calcula a função de onda de modo explícito
        u_1 = Calcula_funcao_onda(alpha, beta, c2, delta_x, nt, T, x_c, delta_t, nx)

        u  = []
        for item in xk:
            u_item = Calcula_funcao_onda(alpha, beta, c2, delta_x, nt, tf, item, delta_t, nx)
            u.append(u_item)

        B = monta_matriz_B(u, ti, tf, xk, xr, delta_x, delta_t, K)
        C = monta_matriz_C(u, ti, tf, xk, xr, delta_x, delta_t, K)
        print("passei??")
        a = [0, 0, 0]
        a = metodo_SOR(B, a, C, w = 1.6, number_it = 10000)
        print(a)
        #Criação do eixo x da função de plotagem
        eixo_x = []
        for j in range(0, nx+1):
            xj = j*delta_x
            eixo_x.append(xj)

        #Plotagem da função de onda em um determinado tempo t
        while True:
            print("Escolha o tempo de plotagem [0~", end='')
            print(T, end='')
            print("]: ", end='')
            tempo = float(input())
            plota_grafico(u_1, tempo, eixo_x, delta_t, nt)
            resp = input("Gostaria de plotar para outro tempo [s/n]: ")
            if resp != "s":
                break
        break

#-------------------------------------------------------------------------------
def F_Força(c2, alpha, beta, t, x, x_c):
    '''(float, float, int, float, float, float,) -> (float)

    A função a seguir recebe alguns parâmetros c2, alpha, beta, t,
    x, x_c e retorna o valor de fonte calcula em x. Caso x != x_c
    a função retorna 0.
    '''

    if (abs(x - x_c) > 0.0001):
        return 0
    termo = (beta**2)*(t**2)*(math.pi**2)
    f_valor = 1000*c2*(1 - 2*termo)*(math.e**(-termo))

    return f_valor

#-------------------------------------------------------------------------------
def Calcula_funcao_onda(alpha, beta, c2, delta_x, nt, T, x_c, delta_t, nx):
    '''(float, float, int, float, int, float, float, float, int) -> (None)

    A função a seguir recebe alguns parâmetros alpha, beta, c2, delta_x,
    nt, T, x_c, delta_t, nx e retorna uma lista u com a função de onda
    a discretizada ao passo de delta_x e delta_t.
    '''

    u =[]
    #Inicia a matriz zerada
    for i in range(nt+1):
        val_col = []
        for j in range(nx+1):
            val_col.append(0)
        u.append(val_col)

    #Calcula os valores de u[i+1][j]
    for i in range(1, nt):
        ti = i*delta_t

        for j in range(0, nx):
            xj = j*delta_x

            #Condições de Dirichlet
            if j == 0 or j == nx:
                u[i][j] = 0
                continue

            #Cálculo do valor da fonte
            fx = F_Força(c2, alpha, beta, ti, xj, x_c)

            #Cálculo de u
            alpha_2 = alpha**2
            delta_t_2 = delta_t**2
            u[i+1][j] = -u[i-1][j] + 2*(1-alpha_2)*u[i][j] + alpha_2*(u[i][j+1] + u[i][j-1]) + delta_t_2*fx

    return u

#-------------------------------------------------------------------------------
def plota_grafico(u, t, eixo_x, delta_t, nt):
    '''(list, float, list, float, int) -> (None)

    A função a seguir recebe uma lista u de valores de uma
    função de onda, um tempo t, um eixo x e um delta_t e
    plota o gráfico de u em um determinado tempo t.
    '''

    eixo_y = u[int(t/delta_t)]
    plt.plot(eixo_x, eixo_y)
    plt.xlabel("x")
    plt.ylabel("u(t,x)")
    title = "Solução equação onda (n_t=" + str(nt) + ", t=" + str(t) + ")"
    plt.title(title)
    plt.ylim(-0.6, 0.6)
    plt.show()

#-------------------------------------------------------------------------------
def monta_matriz_B(u, ti, tf, xk, xr, delta_x, delta_t, K):

    B = []
    for i in range(K):
        B_sub = []
        for j in range (K):
            B_sub.append(-1)
        B.append(B_sub)

    for i in range(K):
        for j in range (K):
            if B[i][j] == -1:
                l = int(xr/delta_x)
                Somatorio = 0
                ini = int(ti/delta_t) + 1
                fim = int(tf/delta_t)
                for k in range(ini, fim):
                    Somatorio += u[i][k][l]*u[j][k][l]
                g0 = u[i][ini-1][l]*u[j][ini-1][l]
                gn = u[i][fim][l]*u[j][fim][l]
                B[i][j] = (delta_t/2)*(g0 + gn + 2*Somatorio)
                B[j][i] = B[i][j]
    return B

#-------------------------------------------------------------------------------
def monta_matriz_C(u, ti, tf, xk, xr, delta_x, delta_t, K):

    C = []
    for i in range(K):
        C.append(0)

    dr3 = numpy.load('dr3.npy')
    for i in range (K):
        l = int(xr/delta_x)
        Somatorio = 0
        ini = int(ti/delta_t) + 1
        fim = int(tf/delta_t)
        for k in range(ini, fim):
            Somatorio += u[i][k][l]*dr3[k]
        g0 = u[i][ini-1][l]*dr3[ini-1]
        gn = u[i][fim][l]*dr3[fim]
        C[i] = (delta_t/2)*(g0 + gn + 2*Somatorio)

    return C

#-------------------------------------------------------------------------------
def metodo_SOR(B, a, C, w, number_it):

    for k in range(number_it):
        for i in range(len(a)):
            Somatorio1 = 0
            j = 0
            while j < i:
                Somatorio1 += B[i][j]*a[j]
                j += 1
            Somatorio2 = 0
            j = i+1
            while j < len(a):
                Somatorio2 += B[i][j]*a[j]
                j += 1
            a[i] = (1-w)*a[i] + (w/B[i][i])*(C[i] - Somatorio1 - Somatorio2)

    return a

#-------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
