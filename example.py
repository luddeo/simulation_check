from simulation_class import MesoRDsimulation

# run each function in the MesoRDsimulation class

new_sim = MesoRDsimulation('simulation_example',0.01)

print("---- Species ----")
print(new_sim.get_species())
new_sim.set_species_order(["Aa","Bb","Cc"])
print(new_sim.get_species())

print()
print("---- Colour ----")
print(new_sim.get_colors())
new_sim.set_colors({"Aa":"#FF0000", "Bb": "#00FF00", "Cc": "#0000FF"})
print(new_sim.get_colors())

print()
print("---- pOcc ----")
print(new_sim.get_pOcc_mean())
print(new_sim.get_pOcc_sd())

print()
print("---- pOcc, threshold = 20 ----")
print(new_sim.get_pOcc_mean(20))
print(new_sim.get_pOcc_sd(20))

print()
print("---- D ----")
print(new_sim.get_D_mean())
print(new_sim.get_D_sd())

print()
print("---- D, threshold = 20 ----")
print(new_sim.get_D_mean(20))
print(new_sim.get_D_sd(20))

print()
print("---- DT ----")
print(new_sim.get_DT_mean())
print(new_sim.get_DT_sd())

print()
print("---- DT, threshold = 20 ----")
print(new_sim.get_DT_mean())
print(new_sim.get_DT_sd())

print()
print("---- plot pOcc ----")
new_sim.plot_pOcc()

print()
print("---- plot trajectory, ID = 30 ----")
new_sim.plot_trajectory(30)

print()
print("---- plot trajectory, ID = 30, lower = 20, upper = 20 ----")
new_sim.plot_trajectory(30, lower_time = 20, upper_time = 30)

print()
print("---- plot trajectory radial, ID = 30 ----")
new_sim.plot_trajectory_radial(30)

print()
print("---- plot trajectory radial, ID = 30, lower = 20, upper = 20 ----")
new_sim.plot_trajectory_radial(30, lower_time = 20, upper_time = 30)

print()
print("---- plot histogram ----")
new_sim.plot_hist()
