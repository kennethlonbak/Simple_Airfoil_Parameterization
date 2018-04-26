import bezier
import pylab as py

def param2coor(x, y, TE_wedge, TE_angle,
               th_tan_LE = 0.5, th_width = 0.5, th_tan_TE = 0.5, th_y_TE = 0.001,
               ca_tan_x = 0.5, ca_tan_LE=0.5,
               coordinate_dir = "counter_clock", N = 101, return_param = False, **kwargs):

    th_CP, th_curve, ca_CP, ca_curve = param2bz_curve(x, y, TE_wedge, TE_angle,
                                                    th_tan_LE, th_width, th_tan_TE, th_y_TE,
                                                    ca_tan_x, ca_tan_LE)


    # Create airfoil from thickness and camber
    x_vec = (1 - py.cos(py.linspace(0, py.pi, N, dtype=float))) / 2
    # Convert bezier to x spaced values
    y_th = bezier2xspaced(th_curve, x_vec)
    y_ca = bezier2xspaced(ca_curve, x_vec)

    # Making upper and lower surface
    y_up = y_th+y_ca
    y_lo = -y_th + y_ca

    x_out = py.concatenate((x_vec[1:][::-1],x_vec))
    if coordinate_dir == "clock":
        y_out = py.concatenate((y_lo[1:][::-1],y_up))
    else:
        y_out = py.concatenate((y_up[1:][::-1],y_lo))

    if return_param:
        import collections
        par_out = collections.OrderedDict()
        names = "x,y,TE_wedge,TE_angle,th_tan_LE,th_width,th_tan_TE,th_y_TE,ca_tan_x,ca_tan_LE".split(",")
        for name, val in locals().items():
            if name in names:
                par_out[name] = val
        return py.array([x_out, y_out]).T, par_out
    else:
        return py.array([x_out,y_out]).T



def param2bz_curve(x, y, TE_wedge, TE_angle,
                   th_tan_LE = 0.5, th_width = 0.5, th_tan_TE = 0.5, th_y_TE = 0.001,
                   ca_tan_LE = 0.8, ca_tan_TE=0.5):

    # Thickness curve
    th_CP = th_param2bz_conp(x, y, TE_wedge,
                             th_tan_LE,
                             th_width,
                             th_tan_TE, th_y_TE)
    th_curve = bezier.Curve(th_CP, degree=len(th_CP[0]) - 1)

    # Camber curve
    ca_CP = ca_param2bz_conp(TE_angle,
                            ca_tan_LE, ca_tan_TE,
                            th_width, x, y)
    ca_curve = bezier.Curve(ca_CP, degree=len(ca_CP[0]) - 1)
    return th_CP, th_curve, ca_CP, ca_curve

def th_param2bz_conp(x, y, TE_wedge,
                     th_tan_LE = 0.5,
                     th_width= 0.1,
                     th_tan_TE = 0.5, th_y_TE = 0.001):
    CP = []
    # LE controle point
    CP.append([0.0,0.0])

    # LE Tangent point
    CP.append([0.0,y*th_tan_LE])

    # Thickness max points
    if x > 0.5:
        th_width_real = (1-x)*th_width
    else:
        th_width_real = x* th_width
    CP.append([x-th_width_real,y])
    CP.append([x+th_width_real,y])

    # TE Tangent point
    CP.append([(x+th_width_real-1)/(y-th_y_TE) *y* th_tan_TE * TE_wedge + 1, y * th_tan_TE+th_y_TE])

    # TE point
    CP.append([1.0,th_y_TE])
    return py.array(CP).T

def ca_param2bz_conp(TE_angle,
                     ca_tan_LE= 0.5, ca_tan_TE=0.5,
                     th_width=0.05, x = 0.5, y = 0.2):
    CP = []
    # LE controle point
    CP.append([0.0, 0.0])

    # LE Tangent point
    CP.append([ca_tan_LE * ca_tan_TE, 0.0])

    # TE Tangent point
    if x > 0.5:
        th_width_real = (1-x)*th_width
    else:
        th_width_real = x* th_width
    CP.append([ca_tan_TE, TE_angle * y / (x + th_width_real - 1) * (ca_tan_TE - 1.0)])

    # TE point
    CP.append([1.0, 0.0])
    return py.array(CP).T

def bezier2xspaced(bz_curve,x_vec):
    s = py.linspace(0,1,len(x_vec)*3)
    coor = bz_curve.evaluate_multi(s)
    y = py.interp(x_vec,coor[0,:],coor[1,:])
    return y

if __name__ == "__main__":
    fig, ax = py.subplots(3,1)
    x = 0.1
    y = 0.1
    TE_wedge = 0.5
    TE_angle = 0.5
    # Create thickness shape
    th_CP, th_curve, ca_CP, ca_curve = param2bz_curve(x, y, TE_wedge, TE_angle)
    th_curve.plot(num_pts=250,ax=ax[0])
    ax[0].plot(*th_CP,'.-')
    ax[0].axis("equal")
    ca_curve.plot(num_pts=250,ax=ax[1])
    ax[1].plot(*ca_CP,'.-')
    ax[1].axis("equal")

    coor, param = param2coor(x, y, TE_wedge, TE_angle, return_param=True)
    print(param)
    ax[2].plot(*coor.T, '-')
    ax[2].axis("equal")
    py.show()



