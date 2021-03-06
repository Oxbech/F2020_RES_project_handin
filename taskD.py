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
hours_in_2015 = pd.date_range("2015-01-01T00:00Z","2015-12-31T23:00Z",
                              freq="H")

# seting the timestamps to calulate for
network.set_snapshots(hours_in_2015)

# adding the electricity system/bus, name: "electricity bus"
network.add("Bus","electricity bus")

# electricity demand in MWh for Austria
demand_elec = data_elec["AUT"]

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
network.add("Carrier","electricity")
network.add("Carrier","hydrogen")


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



## ## ## ## ## ## ## ## CONSTRAIN THE NETWORK ## ## ## ## ## ## ## ##
# ton CO2 equivalents emitted from energy in 1990
co2_1990 = 14e6
# 10% of 1990 emissions allowed
co2_percentage = 0.05
# calculating the equivalent limits on CO2 emissions in ton CO2 equivalents
co2_limit = co2_percentage*co2_1990 # tonCO2e

network.add("GlobalConstraint",
            "co2_limit",
            type="primary_energy",
            carrier_attribute="co2_emissions",
            sense="<=",
            constant=co2_limit
            )


# initializing arrays to store the solutions
p_load = np.zeros([network.snapshots.size,3],dtype="float64")
p_hydro = np.zeros([network.snapshots.size,3],dtype="float64")
p_onshore = np.zeros([network.snapshots.size,3],dtype="float64")
p_solar = np.zeros([network.snapshots.size,3],dtype="float64")
p_OCGT = np.zeros([network.snapshots.size,3],dtype="float64")
p_inverter = np.zeros([network.snapshots.size,3],dtype="float64")
p_electrolysis = np.zeros([network.snapshots.size,3],dtype="float64")
p_fuelCell = np.zeros([network.snapshots.size,3],dtype="float64")


## ## ## ## ## ## ## ## SOLVE THE NETWORK ## ## ## ## ## ## ## ##
# solving the network without any storage
print("Solution 1 of 3: No storage")
network.lopf(network.snapshots,solver_name="gurobi")

mix_noStorage = [network.storage_units_t.p["hydro"].sum(),
                 network.generators_t.p["onshorewind"].sum(),
                 network.generators_t.p["solar"].sum(),
                 network.generators_t.p["OCGT"].sum()]

p_load[:,0] = network.loads_t.p["load"]
p_hydro[:,0] = network.storage_units_t.p["hydro"]
p_onshore[:,0] = network.generators_t.p["onshorewind"]
p_solar[:,0] = network.generators_t.p["solar"]
p_OCGT[:,0] = network.generators_t.p["OCGT"]

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

# co2-emissions in tCO2e
co2_emissions_total = network.generators_t.p["OCGT"].sum() / \
    efficiency_OCGT * network.carriers[:"gas"]["co2_emissions"].mean()
print("\nTotal CO2-emissions: {:1.2f} MtCO2e".format(co2_emissions_total*1e-6))
print("Relative CO2-emissions: {:2.1f}%".format(co2_emissions_total*1e-6 / 
                                                14*100))



## ## ## ## ## ## ## ## ADD BATTERIES ## ## ## ## ## ## ## ##
# annualizing capital (overnight and FOM) costs for battery (inverter) power 
# capacity over 20 years with a discount rate of 0.07
capital_cost_inverter = annuity(20,0.07)*310e3*(1 + 0.01) # €/MW

# capital cost of batteries for every MWh storage capacity, annualized over 
# 15 years with a discount rate of 0.07
capital_cost_batteries = annuity(20,0.07)*144.6e3*(1 + 0) # €/MWh

# battery inverter efficiency
efficiency_inverter = 0.9

# adding a battery bus
network.add("Bus",
            "battery bus",
            carrier="DC")

# adding a link to charge and discharge the batteries from the grid
network.add("Link",
            "inverter",
            bus0="battery bus",
            bus1="electricity bus",
            p_nom_extendable=True,
            p_min_pu=-1, # link is reversible
            capital_cost=capital_cost_inverter*efficiency_inverter,
            marginal_cost=0,
            efficiency=efficiency_inverter
            )

# adding a store to store the electricity in the batteries
network.add("Store",
            "batteries",
            bus="battery bus",
            e_nom_extendable=True,
            capital_cost=capital_cost_batteries,
            e_cyclic=True,
            )

# solving the network with only battery storage
print("\nSolution 2 of 3: Only battery storage")
network.lopf(network.snapshots,solver_name="gurobi")

mix_batteries = [network.storage_units_t.p["hydro"].sum(),
                 network.generators_t.p["onshorewind"].sum(),
                 network.generators_t.p["solar"].sum(),
                 network.generators_t.p["OCGT"].sum()]

p_load[:,1] = network.loads_t.p["load"]
p_hydro[:,1] = network.storage_units_t.p["hydro"]
p_onshore[:,1] = network.generators_t.p["onshorewind"]
p_solar[:,1] = network.generators_t.p["solar"]
p_OCGT[:,1] = network.generators_t.p["OCGT"]
p_inverter[:,1] = network.links_t.p1["inverter"]

print("\nNominal wind power capacity: {:4.0f} MW".format(
    network.generators.p_nom_opt["onshorewind"]))
print("Nominal solar power capacity: {:4.0f} MW".format(
    network.generators.p_nom_opt["solar"]))
print("Nominal OCGT power capacity: {:4.0f} MW".format(
    network.generators.p_nom_opt["OCGT"]))
print("Nominal hydro power capacity: {:4.0f} MW".format(
    network.storage_units.p_nom_opt["hydro"]))

print("\nNominal battery power capacity: {:4.0f} MW".format(
    network.links.p_nom_opt["inverter"]))
print("Nominal battery storage capacity: {:4.0f} MWh".format(
    network.stores.e_nom_opt["batteries"]))

print("\nTotal system cost: {:4.0f} million €".format(
    network.objective*1e-6))
print("Marginal system cost: {:2.1f} €/MWh".format(
    network.objective/network.loads_t.p["load"].sum()))

# co2-emissions in tCO2e
co2_emissions_total = network.generators_t.p["OCGT"].sum() / \
    efficiency_OCGT * network.carriers[:"gas"]["co2_emissions"].mean()
print("\nTotal CO2-emissions: {:1.2f} MtCO2e".format(co2_emissions_total*1e-6))
print("Relative CO2-emissions: {:2.1f}%".format(co2_emissions_total*1e-6 / 
                                                14*100))



## ## ## ## ## ## ## ## ADD HYDROGEN STORAGE ## ## ## ## ## ## ## ##
# annualizing capital (overnight and FOM) costs for hydrogen power capacity,
# e.g. electrolysis and fuel cells, annualized over 18 and 20 years 
# respectively with a discount rate of 0.07
capital_cost_fuelCell = annuity(20,0.07)*339e3*(1 + 0.03) # €/MW

capital_cost_electrolysis = annuity(18,0.07)*350e3*(1 + 0.04) # €/MW

# capital cost of hydrogen storage, for every MWh storage capacity, annualized 
# over 20 years with a discount rate of 0.07
capital_cost_hydrogen = annuity(20,0.07)*8.4e3*(1 + 0) # €/MWh

# fuel cell (hydrogen to AC) efficiency
efficiency_fuelCell = 0.58

# electrolysis (AC to hydrogen) efficiency
efficiency_electrolysis = 0.8

# adding a battery bus
network.add("Bus",
            "hydrogen bus",
            carrier="hydrogen")

# adding a link to act as the fuel cell
network.add("Link",
            "fuel cell",
            bus0="hydrogen bus",
            bus1="electricity bus",
            p_nom_extendable=True,
            capital_cost=capital_cost_fuelCell*efficiency_fuelCell,
            efficiency=efficiency_fuelCell
            )

# adding a link to act as the elctrolysis
network.add("Link",
            "electrolysis",
            bus0="electricity bus",
            bus1="hydrogen bus",
            p_nom_extendable=True,
            capital_cost=capital_cost_electrolysis*efficiency_electrolysis,
            efficiency=efficiency_electrolysis
            )

# adding a store to store the hydrogen
network.add("Store",
            "hydrogen storage",
            bus="hydrogen bus",
            e_nom_extendable=True,
            capital_cost=capital_cost_hydrogen,
            e_cyclic=True
            )

# solving the network with both storages
print("\nSolution 3 of 3: Battery and hydrogen storage")
network.lopf(network.snapshots,solver_name="gurobi")

mix_storage = [network.storage_units_t.p["hydro"].sum(),
               network.generators_t.p["onshorewind"].sum(),
               network.generators_t.p["solar"].sum(),
               network.generators_t.p["OCGT"].sum()]

p_load[:,2] = network.loads_t.p["load"]
p_hydro[:,2] = network.storage_units_t.p["hydro"]
p_onshore[:,2] = network.generators_t.p["onshorewind"]
p_solar[:,2] = network.generators_t.p["solar"]
p_OCGT[:,2] = network.generators_t.p["OCGT"]
p_inverter[:,2] = network.links_t.p1["inverter"]
p_electrolysis[:,2] = network.links_t.p0["electrolysis"]
p_fuelCell[:,2] = network.links_t.p1["fuel cell"]

print("\nNominal wind power capacity: {:4.0f} MW".format(
    network.generators.p_nom_opt["onshorewind"]))
print("Nominal solar power capacity: {:4.0f} MW".format(
    network.generators.p_nom_opt["solar"]))
print("Nominal OCGT power capacity: {:4.0f} MW".format(
    network.generators.p_nom_opt["OCGT"]))
print("Nominal hydro power capacity: {:4.0f} MW".format(
    network.storage_units.p_nom_opt["hydro"]))

print("\nNominal battery power capacity: {:4.0f} MW".format(
    network.links.p_nom_opt["inverter"]))
print("Nominal battery storage capacity: {:4.0f} MWh".format(
    network.stores.e_nom_opt["batteries"]))

print("\nNominal electrolysis power capacity: {:4.0f} MW".format(
    network.links.p_nom_opt["electrolysis"]))
print("Nominal fuel cell power capacity: {:4.0f} MW".format(
    network.links.p_nom_opt["fuel cell"]))
print("Nominal hydrogen storage capacity: {:4.0f} MWh".format(
    network.stores.e_nom_opt["hydrogen storage"]))

print("\nTotal system cost: {:4.0f} million €".format(
    network.objective*1e-6))
print("Marginal system cost: {:2.1f} €/MWh".format(
    network.objective/network.loads_t.p["load"].sum()))

# co2-emissions in tCO2e
co2_emissions_total = network.generators_t.p["OCGT"].sum() / \
    efficiency_OCGT * network.carriers[:"gas"]["co2_emissions"].mean()
print("\nTotal CO2-emissions: {:1.2f} MtCO2e".format(co2_emissions_total*1e-6))
print("Relative CO2-emissions: {:2.1f}%".format(co2_emissions_total*1e-6 / 
                                                14*100))


#%%
## ## ## ## ## ## ## ## DATA-PROCESSING ## ## ## ## ## ## ## ##
# y-limits for dispatch time series
ylim_min = -10000
ylim_max = 20000

# plotting the dispatch time series for three days in summer
startdate = "2015-06-07"
enddate = "2015-06-10"
fig = plt.figure()
plt.plot(network.loads_t.p["load"].loc[startdate:enddate], 
         color="black", linestyle="dashed", label="demand")
plt.plot(network.storage_units_t.p["hydro"].loc[startdate:enddate], 
         color="blue", label="hydro")
plt.plot(network.generators_t.p["onshorewind"].loc[startdate:enddate], 
         color="black", label="onshore wind")
plt.plot(network.generators_t.p["solar"].loc[startdate:enddate], 
         color="orange", label="solar")
plt.plot(network.generators_t.p["OCGT"].loc[startdate:enddate], 
         color="brown", label="gas (OCGT)")
plt.plot(network.links_t.p1["inverter"].loc[startdate:enddate], 
         color="cyan", label="batteries")
plt.plot((-network.links_t.p1["fuel cell"] 
          -network.links_t.p0["electrolysis"]).loc[startdate:enddate], 
         color="purple", label="hydrogen")
fig.autofmt_xdate()
plt.xlim([pd.Timestamp(startdate), pd.Timestamp(enddate)])
plt.ylim([ylim_min,ylim_max])
plt.ylabel("MW electricity")
plt.grid(True)
plt.legend(loc="center", bbox_to_anchor=(.5,1.1), frameon=False, ncol=4)
plt.savefig("../LaTeX/figures/D_disp.eps")


# plotting the annual electricity mix for the three cases
fig = plt.figure()
i = 0; sumHydro = np.asarray([mix_noStorage[i], mix_batteries[i], mix_storage[i]])*1e-6
i = 1; sumWind = np.asarray([mix_noStorage[i], mix_batteries[i], mix_storage[i]])*1e-6
i = 2; sumSolar = np.asarray([mix_noStorage[i], mix_batteries[i], mix_storage[i]])*1e-6
i = 3; sumOCGT = np.asarray([mix_noStorage[i], mix_batteries[i], mix_storage[i]])*1e-6

ind = np.arange(3) # x-locations of the bars
width = 0.35 # width of bars

plt.bar(x=ind, height=sumHydro, width=width, color="blue", label="hydro")
plt.bar(x=ind, height=sumWind, width=width, color="black", label="onshore wind", 
        bottom=sumHydro)
plt.bar(x=ind, height=sumSolar, width=width, color="orange", label="solar", 
        bottom=sumHydro+sumWind)
plt.bar(x=ind, height=sumOCGT, width=width, color="brown", label="gas (OCGT)", 
        bottom=sumHydro+sumWind+sumSolar)
plt.ylabel('TWh electricity')
plt.ylim([0,80])
plt.xticks(ind, ('No storage', 'Batteries only', 'Batt+Hydrogen'))
plt.legend(loc="center", bbox_to_anchor=(.5,1.1), frameon=False , ncol=5)
plt.savefig("../LaTeX/figures/D_mix.eps")


# plotting the dispatch time series for the three cases
start_index=(31+28+31+30+31+6)*24
end_index=start_index+3*24+1
snapshots = pd.date_range("2015-01-01T00:00Z", "2015-12-31T23:00Z", freq="H")
titles = ["No storage", "Batteries only", "Batt.+hydrogen"]
fileNames = ["noStorage", "batteries", "both"]

for i in [0,1,2]:
    fig = plt.figure()
    plt.plot(snapshots[start_index:end_index], 
             p_load[start_index:end_index,i], color='black', label='demand', 
             linestyle="dashed")
    plt.plot(snapshots[start_index:end_index], 
             p_hydro[start_index:end_index,i], color="blue", label="hydro")
    plt.plot(snapshots[start_index:end_index], 
             p_onshore[start_index:end_index,i], color="black", label="onshore wind")
    plt.plot(snapshots[start_index:end_index], 
             p_solar[start_index:end_index,i], color="orange", label="solar")
    plt.plot(snapshots[start_index:end_index], 
             p_OCGT[start_index:end_index,i], color='brown', label='gas (OCGT)')
    plt.plot(snapshots[start_index:end_index], 
             p_inverter[start_index:end_index,i], color="cyan", label="batteries")
    plt.plot(snapshots[start_index:end_index], 
             (-p_fuelCell-p_electrolysis)[start_index:end_index,i], 
             color="purple", label="hydrogen")
    plt.xlim([snapshots[start_index], snapshots[end_index-1]])
    fig.autofmt_xdate()
    #plt.title("{}".format(titles[i]))
    plt.ylim([-10000,25000])
    plt.ylabel("MW electricity")
    plt.grid(True)
    plt.legend(loc="center", bbox_to_anchor=(.5,1.1), frameon=False, ncol=4)
    plt.savefig("../LaTeX/figures/D_disp_{}.eps".format(fileNames[i]),
                bbox_inches = "tight")