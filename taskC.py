# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 14:04:30 2020

@author: Anders
"""

## ## ## ## ## ## ## ## IMPORTING PACKAGES ## ## ## ## ## ## ## ##
import pypsa
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



## ## ## ## ## ## ## ## IMPORTING DATA ## ## ## ## ## ## ## ##
# electricity demand in MWh
data_elec = pd.read_csv("../data/electricity_demand.csv", sep=";", index_col=0)

# onshore wind generation capacity factors
data_wind = pd.read_csv("../data/CF_onshore_wind.csv", sep=";", index_col=0)

# PV capacity factors
data_solar = pd.read_csv("../data/CF_pv_optimal.csv", sep=";", index_col=0)

# Hydro inflow in GWh
data_hydro = pd.read_csv("../data/hydro_inflow.csv", sep=",")

# creating a column of dates, by concatenating the three first columns
data_hydro["date"] = data_hydro["Year"].astype(str) + "-" + \
    + data_hydro["Month"].astype(str) + "-" + data_hydro["Day"].astype(str)
# converts the date column to pandas datetime-type
data_hydro["date"] =  pd.to_datetime(data_hydro["date"])
# removes the extraneous columns
data_hydro = pd.DataFrame({"date": data_hydro["date"], 
                           "inflow": data_hydro["Inflow [GWh]"]})
# adding a "Not A Number" row to get data resampled for the 31st of december 
# of the final year of data, 2012
data_hydro.loc[len(data_hydro.index)] = [pd.to_datetime("2013-01-01T00:00")
                                         ,float("NaN")]
# sets the date column to be the index of the DataFrame
data_hydro = data_hydro.set_index("date")
# divides the inflow with 24 hours to get the inflow per hour, i.e. as power 
# in MW
data_hydro["inflow"] = data_hydro["inflow"]/24*1e3
# resamples the data into 1 hour intervals using the foward fill method, i.e. 
# the last valid observation will be propagated forward
data_hydro = data_hydro.resample("H",convention="end").ffill()
# removes all rows with no data, i.e. the final value, which was required to 
# get values for all hours of the 31st of december 2012
data_hydro = data_hydro.dropna(how="any",axis=0)
# creates a series from the DataFrame for the last non-leap year inflow data 
# is available for, 2011 [MW]
hydro_inflow = data_hydro.loc["2011-01-01":"2011-12-31","inflow"]



## ## ## ## ## ## ## ## FUNCTIONS ## ## ## ## ## ## ## ##
def annuity(n,r):
    """Calculate the annuity factor for an asset with a liftime of n years and 
    discount rate of r, e.g. annuity(20,0.05)*20 = 1.6"""
    
    if r > 0:
        return r/( 1. - 1. / (1. + r)**n )
    else:
        return 1/n



## ## ## ## ## ## ## ## CREATE NETWORK ## ## ## ## ## ## ## ##
# naming the network "network"
network = pypsa.Network()

# creating series of hours to be included in the calculation
hours_in_x = pd.date_range("2015-01-01T00:00Z","2015-12-31T23:00Z",
                              freq="H")

# seting the timestamps to calulate for
network.set_snapshots(hours_in_x)

# adding the electricity system/bus, name: "electricity bus"
network.add("Bus","electricity bus")

# electricity demand in MWh for Austria, reset to make the index match the 
# snapshots included in the network
demand_elec = data_elec["AUT"]
demand_elec = demand_elec.reset_index(drop=True)
demand_elec = pd.DataFrame({"date": network.snapshots, "load": demand_elec})
demand_elec = demand_elec.set_index("date")
demand_elec = demand_elec.loc[:,"load"]

# adding the load to the electricty bus, i.e. the electricity demand. 
# Name: "load"
network.add("Load","load",bus="electricity bus",p_set=demand_elec)



## ## ## ## ## ## ## ## ADD GENERATORS ## ## ## ## ## ## ## ##
# onshore wind generation capacity factors
CF_wind = data_wind["AUT"]\
    [[hour.strftime("%Y-%m-%dT%H:%M:%SZ") for hour in network.snapshots]]

# PV capacity factors
CF_solar = data_solar["AUT"]\
    [[hour.strftime("%Y-%m-%dT%H:%M:%SZ") for hour in network.snapshots]]

# changing the index of the hydro inflow series to match the dates of the 
# network snapshots.
hydro_inflow = hydro_inflow.reset_index(drop=True)
hydro_inflow = pd.DataFrame({"date": network.snapshots, "inflow": hydro_inflow})
hydro_inflow = hydro_inflow.set_index("date")
hydro_inflow = hydro_inflow.loc[:,"inflow"]


# adding the different energy carriers, only gas emits CO2
network.add("Carrier","gas",co2_emissions=0.19) # ton_CO2/MWh_th
network.add("Carrier","onshorewind")
network.add("Carrier","solar")
network.add("Carrier","hydro")


# annualizing capital (overnight and FOM) costs for onshore wind over 30 years 
# with a discount rate of 0.07
capital_cost_onshorewind = annuity(30,0.07)*910e3*(1 + 0.033) # €/MW

# adding onshore wind generator
network.add("Generator",
            "onshorewind",
            bus="electricity bus",
            p_nom_extendable=True, # the capacity can be extended
            carrier="onshorewind",
            capital_cost=capital_cost_onshorewind,
            marginal_cost=0, # no fuel costs
            p_max_pu=CF_wind
            )


# annualizing capital (overnight and FOM) costs for utility scale solar 
# generation over 25 years with a discount rate of 0.07
capital_cost_solar = annuity(25,0.07)*425e3*(1 + 0.03) # €/MW

# adding ultility scale solar power generator
network.add("Generator",
            "solar",
            bus="electricity bus",
            p_nom_extendable=True, # the capacity can be extended
            carrier="solar",
            capital_cost=capital_cost_solar,
            marginal_cost=0, # no fuel costs
            p_max_pu=CF_solar
            )


# annualizing capital (overnight and FOM) costs for OCGT, Open Cycle Gas 
# Turbine which serves as backup-power, over 25 years with a discount rate of 
# 0.07
capital_cost_OCGT = annuity(25,0.07)*560e3*(1 + 0.033) # €/MW

# fuel cost of gas for OCGT in €/MWh_th
fuel_cost_OCGT = 21.6

# efficiency of OCGT
efficiency_OCGT = 0.39

# marginal (fuel) cost for OCGT per MWh electricity produced, €/MWh_el
marginal_cost_OCGT = fuel_cost_OCGT / efficiency_OCGT

# adding OCGT (backup) power generator
network.add("Generator",
            "OCGT",
            bus="electricity bus",
            p_nom_extendable=True, # the capacity can be extended
            carrier="gas",
            capital_cost=capital_cost_OCGT,
            marginal_cost=marginal_cost_OCGT,
            efficiency=efficiency_OCGT
            )


# annualizing capital (overnight and FOM) costs for run-of-river hydro 
# generation over 80 years with a discount rate of 0.07
capital_cost_hydro = annuity(80,0.07)*2000e3*(1 + 0.01) # €/MW

# adding hydro generation as a storage unit
network.add("StorageUnit",
            "hydro",
            bus="electricity bus",
            p_nom=13427.4, # from Wagner et al. (total)
            p_nom_extendable=False,
            carrier="hydro",
            capital_cost=capital_cost_hydro,
            marginal_cost=0, # no fuel costs
            efficiency_store=0.87,
            efficiency_dispatch=0.87,
            inflow=hydro_inflow,
            p_min_pu=-0.5844, # equivalent to the installed power capacity
            max_hours=3.5065, # equivalent to the installed storage capacity
            cyclic_state_of_charge=True # intial SoC = final SoC
            )



## ## ## ## ## ## ## ## SOLVE THE NETWORK ## ## ## ## ## ## ## ##
# making a list of (non-leap) years to calculate for
years = np.arange(1995,2017,2)
#years = np.arange(2015,2017,2)

# initializing arrays to store the solutions
p_hydro_avg = np.empty([years.size,1],dtype="float64")
p_onshore_avg = np.empty([years.size,1],dtype="float64")
p_solar_avg = np.empty([years.size,1],dtype="float64")
p_OCGT_avg = np.empty([years.size,1],dtype="float64")

# calculating for the various years
for i in range(years.size):
    print("\nSolution {} of {}".format(i+1,years.size))
    
    # creating series of hours to be included in the calculation
    hours_in_x = pd.date_range("{:d}-01-01T00:00Z".format(years[i]),
                               "{:d}-12-31T23:00Z".format(years[i]),
                               freq="H")
    # seting the timestamps to calulate for
    network.set_snapshots(hours_in_x)
    
    # changing the index of the electricity demand series to match the dates of 
    # the network snapshots.
    demand_elec = data_elec["AUT"]
    demand_elec = demand_elec.reset_index(drop=True)
    demand_elec = pd.DataFrame({"date": network.snapshots, "load": demand_elec})
    demand_elec = demand_elec.set_index("date")
    network.loads_t.p_set = demand_elec
    
    # changing the index of the hydro inflow series to match the dates of the 
    # network snapshots.
    hydro_inflow = data_hydro.loc["2011-01-01":"2011-12-31","inflow"]
    hydro_inflow = hydro_inflow.reset_index(drop=True)
    hydro_inflow = pd.DataFrame({"date": network.snapshots, "hydro": hydro_inflow})
    hydro_inflow = hydro_inflow.set_index("date")
    network.storage_units_t.inflow = hydro_inflow
    
    # onshore wind generation capacity factors for the current year
    CF_wind = data_wind["AUT"]\
        [[hour.strftime("%Y-%m-%dT%H:%M:%SZ") for hour in network.snapshots]]
    network.generators_t.p_max_pu["onshorewind"] = CF_wind
    
    # PV capacity factors for the current year
    CF_solar = data_solar["AUT"]\
        [[hour.strftime("%Y-%m-%dT%H:%M:%SZ") for hour in network.snapshots]]
    network.generators_t.p_max_pu["solar"] = CF_solar
    
    # solving the network for the current year
    network.lopf(network.snapshots,solver_name="gurobi")
    
    # storing the solution for the current year
    p_hydro_avg[i] = network.storage_units_t.p["hydro"].mean()
    p_onshore_avg[i] = network.generators_t.p["onshorewind"].mean()
    p_solar_avg[i] = network.generators_t.p["solar"].mean()
    p_OCGT_avg[i] = network.generators_t.p["OCGT"].mean()
    
    print("\nNominal wind power capacity: {:4.0f} MW".format(
        network.generators.p_nom_opt["onshorewind"]))
    print("Nominal solar power capacity: {:4.0f} MW".format(
        network.generators.p_nom_opt["solar"]))
    print("Nominal OCGT power capacity: {:4.0f} MW".format(
        network.generators.p_nom_opt["OCGT"]))
    print("Nominal hydro power capacity: {:4.0f} MW".format(
        network.storage_units.p_nom_opt["hydro"]))
    
    print("\nTotal system cost: {:4.0f} million €".format(
        network.objective*1e-6))
    print("Marginal system cost: {:2.1f} €/MWh".format(
        network.objective/network.loads_t.p["load"].sum()))
    

#%%
## ## ## ## ## ## ## ## DATA-PROCESSING ## ## ## ## ## ## ## ##
fig = plt.figure()
plt.subplot(221)
plt.plot(years, p_hydro_avg, marker="o", linestyle="None", 
         color="blue", label="hydro")
plt.yticks([0, 1000, 2000, 3000])
plt.grid(True)
plt.xticks(np.arange(years[0], years[-1]+1, step=5))
fig.autofmt_xdate()

plt.subplot(222)
plt.plot(years, p_onshore_avg, marker="o", linestyle="None", 
         color="black", label="onshore wind")
plt.yticks([0, 1000, 2000, 3000],[])
plt.grid(True)
plt.xticks(np.arange(years[0], years[-1]+1, step=5))
fig.autofmt_xdate()

plt.subplot(223)
plt.plot(years, p_solar_avg, marker="o", linestyle="None", 
         color="orange", label="solar")
plt.yticks([0, 1000, 2000, 3000])
plt.grid(True)
plt.xticks(np.arange(years[0], years[-1]+1, step=5))
fig.autofmt_xdate()

plt.subplot(224)
plt.plot(years, p_OCGT_avg, marker="o", linestyle="None", 
         color="brown", label="gas (OCGT)")
plt.yticks([0, 1000, 2000, 3000],[])
plt.grid(True)
plt.xticks(np.arange(years[0], years[-1]+1, step=5))
fig.autofmt_xdate()

fig.legend(loc="center", bbox_to_anchor=(.5,0.95), frameon=False , ncol=4)
fig.text(0.02, 0.5, "MW electricity", va="center", rotation="vertical")
#fig.tight_layout()
fig.savefig("../LaTeX/figures/C_avgP.eps",bbox_inches = "tight")

print("\nHydro power capacity variability:  {:5.0f} MW".format(
        p_hydro_avg.std()))
print("Wind power capacity variability:   {:5.0f} MW".format(
        p_onshore_avg.std()))
print("Solar power capacity variability:  {:5.0f} MW".format(
        p_solar_avg.std()))
print("OCGT power capacity variability:   {:5.0f} MW".format(
        p_OCGT_avg.std()))

