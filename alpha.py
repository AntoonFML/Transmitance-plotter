import numpy as np
import matplotlib
import sympy as sp
import control as ct

print("Imports complete.")
print("NumPy version:", np.__version__)
print("Matplotlib version:", matplotlib.__version__)

is_stable = True
next_element = [(2,0), (2,1), (3,0), (4,0)]
unstable_roots = 0
boundary_roots = 0
is_negative = False
x = sp.symbols('x')
row_coeffs = {
    0:x**4,
    1:x**3,
    2:x**2,
    3:x**1,
    4:x**0
}

# input coefficients

tran_L = float(input("b0: "))
print(tran_L)
a3  = float(input("a3: "))
a2 = float(input("a2: "))
a1 = float(input("a1: "))
a0 = float(input("a0: "))
print("s**4 +", a3, "* s**3+", a2, "* s**2+", a1, "* s+", a0)
tran_M = [1, a3, a2, a1, a0]

#define the function to calculate the next element of the table

#limit_result = sp.limit(f, x, 0, dir='+')

def calc_next_element(tab, element):
    if element[1] != 1:
        rows = [element[0]-2, element[0]-1]
        cols = [element[1], element[1]+1]
    else:
        rows = [element[0]-2, element[0]-1]
        cols = [element[1]-1, element[1]+1]

    det_matrix = tab[np.ix_(rows, cols)]
    if det_matrix[1][0] == 0 and det_matrix[1][1] != 0:
        det_matrix[1][0] = x
        det_function = (det_matrix[0][0]*det_matrix[1][1]) - (det_matrix[0][1]*det_matrix[1][0])
        det = det / (-det_matrix[1][0])
        det = sp.limit(det_function, x, 0, dir='+')


    # TODO: check for zero row , derivative of the polynomial of the row above, coefficients in place of the zeros, roots of the polynomial
    elif det_matrix[1][0] == 0 and det_matrix[1][1] == 0 and det_matrix[0][1] != 0:
        if element[0] == 0:
            poly = element[0][0] * row_coeffs.get(element[0]) + element[0][1] * (row_coeffs.get(element[0]+2)) + element[0][2] * (row_coeffs.get(element[0]+4))
        else:
            poly = element[0][0] * row_coeffs.get(element[0]) + element[0][1] * (row_coeffs.get(element[0]+2))
        deriv = sp.diff(poly, x)
        sp.Poly(deriv, x)
        coeffs = sp.Poly(deriv, x).all_coeffs()
        coeffs = [float(i) for i in coeffs]
        det_matrix[1][0] = coeffs[0]
        det_matrix[1][1] = coeffs[1]
        det = np.linalg.det(det_matrix)
        det = det / (-det_matrix[1][0])
        roots = sp.all_roots(deriv, x) 
        for i in roots:
            if np.real(i) > 0:
                unstable_roots += 1
                is_stable = False
            elif np.real(i) == 0 and np.imag(i) != 0:
                is_stable = False
                boundary_roots += 1
    else:
        det = np.linalg.det(det_matrix)
        det = det / (-det_matrix[1][0])
    return det

#check for negative coefficients in the characteristic polynomial

for nums in tran_M:
    if nums < 0:
        is_stable = False
        print("The system is unstable, cannot simulate the output.")
        break

det1, det2, det3, det4 = 0, 0, 0, 0
dets = [det1, det2, det3, det4]

tab_routh = np.array([[1, a2, a0],
                     [a3,a1,0],
                     [dets[0], dets[1], 0],
                     [dets[2], 0, 0],
                     [dets[3], 0, 0]])

for index in range(len(next_element)):
    dets[index] = calc_next_element(tab_routh, next_element[index])
    tab_routh = np.array([[1, a2, a0],
                          [a3,a1,0],
                          [dets[0], dets[1], 0],
                          [dets[2], 0, 0],
                          [dets[3], 0, 0]])

print(tab_routh)

first_col = tab_routh[:, 0]
print("First column:", first_col)
for nums in range(1, len(first_col)):
    if first_col[nums]*first_col[nums-1] < 0:
        is_stable = False
        unstable_roots += 1
print("Unstable roots:", unstable_roots)
print("Boundary roots:", boundary_roots)

#the result of the Routh-Hurwitz test
if is_stable:
    print("The system is stable.")
elif unstable_roots == 0 and boundary_roots != 0:
    print("The system is on the boundary of stability.")
else:   
    print("The system is unstable.")

    #TODO: bode plots of the function, graph of input (sine and unit step), graph of output 