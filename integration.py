from deflection_y_axis import x, force_span, tq_span
from spacing import vspace_list, x_list, spacing_x
import numpy as np
from matplotlib import pyplot as plt

def int_w(value):
    coeff_a = []
    coeff_b = []

    # -------------- Find interpolation coefficients ----------

    for i in range(0, 40):
        a = (force_span[i + 1] - force_span[i]) / spacing_x[i]
        b = force_span[i] - a * x[i]
        coeff_a.append(a)
        coeff_b.append(b)

    # -------------- Let's make the mesh finer -----------------

    mesh_fine = []
    mesh = 50
    x_fine = []
    inc_list = []

    for j in range(0, 40):

        mesh_fine.append(force_span[j])  # Append first force

        inc = spacing_x[j] / mesh  # Calculate the increment
        inc_list.append(inc)

        x_fine.append(x[j])

        for w in range(1, mesh + 1):  # Add another element
            mesh_fine.append(coeff_a[j] * (x[j] + inc * w) + coeff_b[j])
            x_fine.append(x[j] + inc * w)

    # ---------- To verify if the interpolation was successful, graphs should overlap ----------

    # plt.plot(x, force_span, label = "Original")
    # plt.plot(x_fine, mesh_fine, label="Fine")
    # plt.legend()
    # plt.show()


    # --------- Start integration for w ---------
    nr_int = value
    integration = [0]

    for j in range(0, nr_int):
        total = 0
        storage = [0]
        if j == 0:
            for i in range(0, len(mesh_fine) - 1):
                total += (x_fine[i + 1] - x_fine[i]) * (mesh_fine[i + 1] + mesh_fine[i]) * 0.5
                integration.append(total)
            # plt.plot(x_fine, mesh_fine, label = "original")
            # plt.pause(0.05)

        else:
            for i in range(0, len(integration) - 1):
                total += (x_fine[i + 1] - x_fine[i]) * (integration[i + 1] + integration[i]) * 0.5
                storage.append(total)
            integration = storage
        name = "Integration nr" + str(j+1)
    #     plt.plot(x_fine, integration, label = name)
    #     plt.pause(0.05)

    # plt.legend()
    # plt.show()

    return integration, x_fine


def int_tau(value):
    # -------- Start integration for tau --------

    coeff_a = []
    coeff_b = []

    # -------------- Find interpolation coefficients ----------

    for i in range(0, 40):
        a = (tq_span[i + 1] - tq_span[i]) / spacing_x[i]
        b = tq_span[i] - a * x[i]
        coeff_a.append(a)
        coeff_b.append(b)

    # -------------- Let's make the mesh finer -----------------

    mesh_fine = []
    mesh = 50
    x_fine = []
    inc_list = []

    for j in range(0, 40):

        mesh_fine.append(tq_span[j])  # Append first force

        inc = spacing_x[j] / mesh  # Calculate the increment
        inc_list.append(inc)

        x_fine.append(x[j])

        for w in range(1, mesh + 1):  # Add another element
            mesh_fine.append(coeff_a[j] * (x[j] + inc * w) + coeff_b[j])
            x_fine.append(x[j] + inc * w)

    # ---------- To verify if the interpolation was successful, graphs should overlap ----------

    # plt.plot(x, force_span, label = "Original")
    # plt.plot(x_fine, mesh_fine, label="Fine")
    # plt.legend()
    # plt.show()

    # --------- Start integration for w ---------
    nr_int = value
    integration = [0]

    for j in range(0, nr_int):
        total = 0
        storage = [0]
        if j == 0:
            for i in range(0, len(mesh_fine) - 1):
                total += (x_fine[i + 1] - x_fine[i]) * (mesh_fine[i + 1] + mesh_fine[i]) * 0.5
                integration.append(total)
            plt.plot(x_fine, mesh_fine, label="original")
            plt.pause(0.05)

        else:
            for i in range(0, len(integration) - 1):
                total += (x_fine[i + 1] - x_fine[i]) * (integration[i + 1] + integration[i]) * 0.5
                storage.append(total)
            integration = storage
        name = "Integration nr" + str(j + 1)
        plt.plot(x_fine, integration, label=name)
        plt.pause(0.05)

    plt.legend()
    plt.show()

    print("The integration for tau is ", integration[-1])