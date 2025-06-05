# Antoni Rumowski 193087 Marek Marcinko 197870

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import control as ct

w_max = 100 # maximum frequency for the Bode plot
dw = 0.001 # frequency step for the Bode plot
N = 4 #order of the system
h = 0.001 # time step for the simulation
T = 10 # total time for the simulation
L = 2.5 # number of periods for the sine response
M = 8 #amplitude of the sine response




print("Imports complete.")
print("NumPy version:", np.__version__)
print("Matplotlib version:", matplotlib.__version__)
print("Control version:", ct.__version__)
input("Press Enter to continue...")

is_stable = True

# input coefficients

tran_numerator = float(input("b0: "))
a3  = float(input("a3: "))
a2 = float(input("a2: "))
a1 = float(input("a1: "))
a0 = float(input("a0: "))
print("s**4 +", a3, "* s**3+", a2, "* s**2+", a1, "* s+", a0)
tran_denominator = [1, a3, a2, a1, a0]

def bode_plot(numerator, denominator):
    w = np.arange(0, w_max, dw)  # frequency range for the Bode plot
    denominator_jw = (complex(0,1)*w)**N + a3*(complex(0,1)*w)**(N-1) + a2*(complex(0,1)*w)**(N-2) + a1*(complex(0,1)*w) + a0
    h_jw = numerator / denominator_jw  # transfer function H(jw)

    Aw = np.abs(h_jw)  # magnitude
    Fw = np.angle(h_jw)  # phase
    for k in range(len(Fw)):
        if Fw[k] < 0:
            Fw[k] += 2 * np.pi

    plt.figure()
    plt.semilogx(w, Aw, 'r-', label='Magnitude')
    plt.grid(True)
    plt.semilogx(w, Fw, 'b-', label='Phase')
    plt.title('Bode Plot')
    plt.xlabel('Frequency (rad/s)')
    plt.legend()
    plt.savefig("bode_plot_together.png")
    plt.close()

    plt.figure()
    plt.subplot(2, 1, 1)
    plt.semilogx(w, 20 * np.log10(Aw))
    plt.title('Bode Plot')
    plt.ylabel('Magnitude (dB)')
    plt.grid(True)

    plt.subplot(2, 1, 2)
    plt.semilogx(w, np.degrees(Fw))
    plt.xlabel('Frequency (rad/s)')
    plt.ylabel('Phase (deg)')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("bode_plot.png")
    plt.close()
    print("Bode plot complete.")

def ct_bode_plot(G):

    mag, phase, omega = ct.bode(G, plot=False)
    for i in range(len(phase)):
        if phase[i] < 0:
            phase[i] += 2 * np.pi
    plt.figure()
    plt.semilogx(omega, mag, 'r-', label='Magnitude')
    plt.grid(True)
    plt.semilogx(omega, phase, 'b-', label='Phase')
    plt.title('Bode Plot')
    plt.xlabel('Frequency (rad/s)')
    plt.savefig("ct_bode_plot_together.png")
    plt.close()

    plt.figure()
    plt.subplot(2, 1, 1)
    plt.semilogx(omega, 20 * np.log10(mag)) 
    plt.title('Bode Plot')
    plt.ylabel('Magnitude (dB)')
    plt.grid(True)

    plt.subplot(2, 1, 2)
    plt.semilogx(omega, np.degrees(phase))
    plt.xlabel('Frequency (rad/s)')
    plt.ylabel('Phase (deg)')
    plt.grid(True)

    plt.tight_layout()
    plt.savefig("ct_bode_plot.png")
    plt.close()
    

def step_response(numerator, denominator):
    # time vector for the step response
    t = np.arange(0, T, h)
    # step input
    u = np.ones_like(t)

    b = [numerator]  # numerator coefficients
    a = denominator  # denominator coefficients

    y = np.zeros_like(t) # initialize output 

    A_matrix = np.matrix([[0, 1, 0, 0],
                          [0, 0, 1, 0],
                          [0, 0, 0, 1],
                          [-a[4], -a[3], -a[2], -a[1]]])
    B_matrix = np.matrix([[0],
                          [0],
                          [0],
                          [1]])
    C_matrix = np.matrix([[b[0], 0, 0, 0]])
    D_matrix = np.matrix([[0]])


    total = np.size(u)/np.size(u[0])
    xi_1 =np.matrix([[0],
                     [0],
                     [0],
                     [0]])

    for i in range(int(total)):
        Ax = A_matrix * xi_1
        Bu = B_matrix * u[i]
        Cx = C_matrix * xi_1
        Du = D_matrix * u[i]
        xi = Ax + Bu
        xi = xi * h
        xi = xi_1 + xi
        xi_1 = xi
        y[i] = (Cx + Du)[0,0]
        

    plt.figure()
    plt.plot(t, y, label="y(t) – step response")
    plt.plot(t, u, '--', label="u(t) – unit step")
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude")
    plt.grid()
    plt.title("Step response of the system")
    plt.legend()
    plt.savefig("step_response.png")
    plt.close()
    print("Step response complete.")

def sine_response(numerator, denominator):
    # time vector for the step response
    t = np.arange(0, T, h)
    u = np.zeros_like(t)  # initialize input
    w = 2 * np.pi * L / T  # sine frequency based on the number of periods
    

    b = [numerator]  # numerator coefficients
    a = denominator  # denominator coefficients

    y = np.zeros_like(t) # initialize output 

    A_matrix = np.matrix([[0, 1, 0, 0],
                          [0, 0, 1, 0],
                          [0, 0, 0, 1],
                          [-a[4], -a[3], -a[2], -a[1]]])
    B_matrix = np.matrix([[0],
                          [0],
                          [0],
                          [1]])
    C_matrix = np.matrix([[b[0], 0, 0, 0]])
    D_matrix = np.matrix([[0]])

  
    xi_1 =np.matrix([[0],
                     [0],
                     [0],
                     [0]])
    for i in range(T/h):
        us = M * np.sin(w * i * h)
        u[i] = us 
        Ax = A_matrix * xi_1
        Bu = B_matrix * us
        Cx = C_matrix * xi_1
        Du = D_matrix * us
        xi = Ax + Bu
        xi = xi * h
        xi = xi_1 + xi
        xi_1 = xi
        y[i] = (Cx + Du)[0,0]
        

    plt.figure()
    plt.plot(t, y, label="y(t) – Sine response")
    plt.plot(t, u, '--', label="u(t) – sine input")
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude")
    plt.grid()
    plt.title("Sine response of the system")
    plt.legend()
    plt.savefig("sine_response.png")
    plt.close()
    print("Sine response complete.")


def ct_sine_response(G):
    t = np.arange(0, T, h)
    u = np.zeros_like(t)  # initialize input
    w = 2 * np.pi * L / T  # sine frequency based on the number of periods
    for i in range(T/h):
        us = M * np.sin(w * i * h)
        u[i] = us



    t_out, y_out = ct.forced_response(G, T=t, U=u)

    plt.figure()
    plt.plot(t_out, u, '--', label="u(t) – sine input")
    plt.plot(t_out, y_out, label="y(t) – system response")
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude")
    plt.legend()
    plt.grid()
    plt.title("System response to sine input")
    plt.savefig("ct_sine_response.png")
    plt.close()
    print("Sine response complete.")

def check_stability(tab):
    rows2 = [0,1]
    cols2 = [0,1]
    matrix2 = tab[np.ix_(rows2, cols2)]

    test2 = np.linalg.det(matrix2)

    rows3 = [0,1,2]
    cols3 = [0,1,2]
    matrix3 = tab[np.ix_(rows3, cols3)]

    test3 = np.linalg.det(matrix3)
    test4 = np.linalg.det(tab)

    print("Matrix 2x2:", matrix2)
    print("Matrix 3x3:", matrix3)
    if test2 <= 0 or test3 <= 0 or test4 <= 0:
        return False
    else:
        return True
    
#check for negative coefficients in the characteristic polynomial

for nums in tran_denominator:
    if nums < 0:
        is_stable = False
        print("The system is unstable, cannot simulate the output.")
        break

tab_hurwitz = np.array([[a3,a1,0,0],
                     [1,a2,a0,0],
                     [0,a3,a1,0],
                     [0,1,a2,a0]])
print(tab_hurwitz)
is_stable = check_stability(tab_hurwitz)

G = ct.TransferFunction([tran_numerator], tran_denominator)

#the result of the Hurwitz test
if is_stable:
    print("The system is stable.")
else:   
    print("The system is unstable.")

#ct_bode_plot(G)
bode_plot(tran_numerator, tran_denominator)
#ct.step_response(G)
step_response(tran_numerator, tran_denominator)
#ct_sine_response(G)
sine_response(tran_numerator, tran_denominator)
input("All tasks complete. Press Enter to exit...")
SystemExit(0)


