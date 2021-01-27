# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 14:04:30 2020

@author: Anders
"""

## ## ## ## ## ## ## ## IMPORTING PACKAGES ## ## ## ## ## ## ## ##
import pypsa
import pandas as pd
import matplotlib.pyplot as plt



## ## ## ## ## ## ## ## IMPORTING DATA ## ## ## ## ## ## ## ##
# electricity demand in MWh
data_elec = pd.read_csv("../data/electricity_demand.csv", sep=";", index_col=0)

# onshore wind generation capacity factors
data_wind = pd.read_csv("../data/CF_onshore_wind.csv", sep=";", index_col=0)

# PV capacity factors
data_solar = pd.read_csv("../data/CF_pv_optimal.csv", sep=";", index_col=0)

# Hydro inflow in GWh
data_hydro_AUT = pd.read_csv("../data/hydro_inflow.csv", sep=",")
data_hydro_SVK = pd.read_csv("../data/hydro_inflow_SVK.csv", sep=",")
data_hydro_CZE = pd.read_csv("../data/hydro_inflow_CZE.csv", sep=",")
data_hydro = pd.DataFrame({"Year": data_hydro_AUT["Year"],
                           "Month": data_hydro_AUT["Month"],
                           "Day": data_hydro_AUT["Day"],
                           "AUT": data_hydro_AUT["Inflow [GWh]"],
                           "SVK": data_hydro_SVK["Inflow [GWh]"],
                           "CZE": data_hydro_CZE["Inflow [GWh]"],
                           })

# creating a column of dates, by concatenating the three first columns
data_hydro["date"] = data_hydro["Year"].astype(str) + "-" + \
    + data_hydro["Month"].astype(str) + "-" + data_hydro["Day"].astype(str)
# converts the date column to pandas datetime-type
data_hydro["date"] =  pd.to_datetime(data_hydro["date"])
# removes the extraneous columns
data_hydro = pd.DataFrame({"date": data_hydro["date"], 
                           "AUT": data_hydro["AUT"],
                           "SVK": data_hydro["SVK"],
                           "CZE": data_hydro["CZE"]
                           })
# adding a "Not A Number" row to get data resampled for the 31st of december 
# of the final year of data, 2012
data_hydro.loc[len(data_hydro.index)] = [pd.to_datetime("2013-01-01T00:00"),
                                         float("NaN"),
                                         float("NaN"),
                                         float("NaN")]
# sets the date column to be the index of the DataFrame
data_hydro = data_hydro.set_index("date")
# divides the inflow with 24 hours to get the inflow per hour, i.e. as power 
# in MW
data_hydro["AUT"] = data_hydro["AUT"]/24*1e3
data_hydro["SVK"] = data_hydro["SVK"]/24*1e3
data_hydro["CZE"] = data_hydro["CZE"]/24*1e3
# resamples the data into 1 hour intervals using the foward fill method, i.e. 
# the last valid observation will be propagated forward
data_hydro = data_hydro.resample("H",convention="end").ffill()
# removes all rows with no data, i.e. the final value, which was required to 
# get values for all hours of the 31st of december 2012
data_hydro = data_hydro.dropna(how="any",axis=0)
# creates a series from the DataFrame for the last non-leap year inflow data 
# is available for, 2011 [MW]
hydro_inflow_AUT = data_hydro.loc["2011-01-01":"2011-12-31","AUT"]
hydro_inflow_SVK = data_hydro.loc["2011-01-01":"2011-12-31","SVK"]
hydro_inflow_CZE = data_hydro.loc["2011-01-01":"2011-12-31","CZE"]



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


# adding the electricity bus for Austria
network.add("Bus","Austrian electricity bus")

# electricity demand in MWh for Austria
demand_AUT = data_elec["AUT"]

# adding the load to the Austrian electricty bus
network.add("Load","Austrian load",bus="Austrian electricity bus",p_set=demand_AUT)


# adding the electricity bus for Slovakia
network.add("Bus","Slovak electricity bus")

# electricity demand in MWh for Slovakia
demand_SVK = data_elec["SVK"]

# adding the load to the Slovak electricty bus
network.add("Load","Slovak load",bus="Slovak electricity bus",p_set=demand_SVK)


# adding the electricity bus for Czechia
network.add("Bus","Czech electricity bus")

# electricity demand in MWh for Czechia
demand_CZE = data_elec["CZE"]

# adding the load to the Czech electricty bus
network.add("Load","Czech load",bus="Czech electricity bus",p_set=demand_CZE)



## ## ## ## ## ## ## ## ADD GENERATORS ## ## ## ## ## ## ## ##
# onshore wind generation capacity factors
CF_wind_AUT = data_wind["AUT"]\
    [[hour.strftime("%Y-%m-%dT%H:%M:%SZ") for hour in network.snapshots]]
CF_wind_SVK = data_wind["SVK"]\
    [[hour.strftime("%Y-%m-%dT%H:%M:%SZ") for hour in network.snapshots]]
CF_wind_CZE = data_wind["CZE"]\
    [[hour.strftime("%Y-%m-%dT%H:%M:%SZ") for hour in network.snapshots]]

# PV capacity factors
CF_solar_AUT = data_solar["AUT"]\
    [[hour.strftime("%Y-%m-%dT%H:%M:%SZ") for hour in network.snapshots]]
CF_solar_SVK = data_solar["SVK"]\
    [[hour.strftime("%Y-%m-%dT%H:%M:%SZ") for hour in network.snapshots]]
CF_solar_CZE = data_solar["CZE"]\
    [[hour.strftime("%Y-%m-%dT%H:%M:%SZ") for hour in network.snapshots]]


# changing the index of the hydro inflow series to match the dates of the 
# network snapshots.
hydro_inflow_AUT = hydro_inflow_AUT.reset_index(drop=True)
hydro_inflow_AUT = pd.DataFrame({"date": network.snapshots, "inflow": hydro_inflow_AUT})
hydro_inflow_AUT = hydro_inflow_AUT.set_index("date")
hydro_inflow_AUT = hydro_inflow_AUT.loc[:,"inflow"]

hydro_inflow_SVK = hydro_inflow_SVK.reset_index(drop=True)
hydro_inflow_SVK = pd.DataFrame({"date": network.snapshots, "inflow": hydro_inflow_SVK})
hydro_inflow_SVK = hydro_inflow_SVK.set_index("date")
hydro_inflow_SVK = hydro_inflow_SVK.loc[:,"inflow"]

hydro_inflow_CZE = hydro_inflow_CZE.reset_index(drop=True)
hydro_inflow_CZE = pd.DataFrame({"date": network.snapshots, "inflow": hydro_inflow_CZE})
hydro_inflow_CZE = hydro_inflow_CZE.set_index("date")
hydro_inflow_CZE = hydro_inflow_CZE.loc[:,"inflow"]


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
            "Austrian onshorewind",
            bus="Austrian electricity bus",
            p_nom_extendable=True, # the capacity can be extended
            carrier="onshorewind",
            capital_cost=capital_cost_onshorewind,
            marginal_cost=0, # no fuel costs
            p_max_pu=CF_wind_AUT
            )
network.add("Generator",
            "Slovak onshorewind",
            bus="Slovak electricity bus",
            p_nom_extendable=True, # the capacity can be extended
            carrier="onshorewind",
            capital_cost=capital_cost_onshorewind,
            marginal_cost=0, # no fuel costs
            p_max_pu=CF_wind_SVK
            )
network.add("Generator",
            "Czech onshorewind",
            bus="Czech electricity bus",
            p_nom_extendable=True, # the capacity can be extended
            carrier="onshorewind",
            capital_cost=capital_cost_onshorewind,
            marginal_cost=0, # no fuel costs
            p_max_pu=CF_wind_CZE
            )


# annualizing capital (overnight and FOM) costs for utility scale solar 
# generation over 25 years with a discount rate of 0.07
capital_cost_solar = annuity(25,0.07)*425e3*(1 + 0.03) # €/MW

# adding ultility scale solar power generator
network.add("Generator",
            "Austrian solar",
            bus="Austrian electricity bus",
            p_nom_extendable=True, # the capacity can be extended
            carrier="solar",
            capital_cost=capital_cost_solar,
            marginal_cost=0, # no fuel costs
            p_max_pu=CF_solar_AUT
            )
network.add("Generator",
            "Slovak solar",
            bus="Slovak electricity bus",
            p_nom_extendable=True, # the capacity can be extended
            carrier="solar",
            capital_cost=capital_cost_solar,
            marginal_cost=0, # no fuel costs
            p_max_pu=CF_solar_SVK
            )
network.add("Generator",
            "Czech solar",
            bus="Czech electricity bus",
            p_nom_extendable=True, # the capacity can be extended
            carrier="solar",
            capital_cost=capital_cost_solar,
            marginal_cost=0, # no fuel costs
            p_max_pu=CF_solar_CZE
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
            "Austrian OCGT",
            bus="Austrian electricity bus",
            p_nom_extendable=True, # the capacity can be extended
            carrier="gas",
            capital_cost=capital_cost_OCGT,
            marginal_cost=marginal_cost_OCGT,
            efficiency=efficiency_OCGT
            )
network.add("Generator",
            "Slovak OCGT",
            bus="Slovak electricity bus",
            p_nom_extendable=True, # the capacity can be extended
            carrier="gas",
            capital_cost=capital_cost_OCGT,
            marginal_cost=marginal_cost_OCGT,
            efficiency=efficiency_OCGT
            )
network.add("Generator",
            "Czech OCGT",
            bus="Czech electricity bus",
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
            "Austrian hydro",
            bus="Austrian electricity bus",
            p_nom=13427.4, # from Wagner et al. (total)
            p_nom_extendable=False,
            carrier="hydro",
            capital_cost=capital_cost_hydro,
            marginal_cost=0, # no fuel costs
            efficiency_store=0.87,
            efficiency_dispatch=0.87,
            inflow=hydro_inflow_AUT,
            p_min_pu=-0.5844, # equivalent to the installed power capacity
            max_hours=3.51, # equivalent to the installed storage capacity
            cyclic_state_of_charge=True # intial SoC = final SoC
            )
network.add("StorageUnit",
            "Slovak hydro",
            bus="Slovak electricity bus",
            p_nom=2031,
            p_nom_extendable=False,
            carrier="hydro",
            capital_cost=capital_cost_hydro,
            marginal_cost=0, # no fuel costs
            efficiency_store=0.87,
            efficiency_dispatch=0.87,
            inflow=hydro_inflow_SVK,
            p_min_pu=-0.4510, # equivalent to the installed power capacity
            max_hours=2.71, # equivalent to the installed storage capacity
            cyclic_state_of_charge=True # intial SoC = final SoC
            )
network.add("StorageUnit",
            "Czech hydro",
            bus="Czech electricity bus",
            p_nom=2260,
            p_nom_extendable=False,
            carrier="hydro",
            capital_cost=capital_cost_hydro,
            marginal_cost=0, # no fuel costs
            efficiency_store=0.87,
            efficiency_dispatch=0.87,
            inflow=hydro_inflow_CZE,
            p_min_pu=-0.5186, # equivalent to the installed power capacity
            max_hours=3.11, # equivalent to the installed storage capacity
            cyclic_state_of_charge=True # intial SoC = final SoC
            )



## ## ## ## ## ## ## ## CONSTRAIN THE NETWORK ## ## ## ## ## ## ## ##
# ton CO2 equivalents emitted from energy in 1990
co2_1990 = 14e6 + 13e6 + 64e6
# 5% of 1990 emissions allowed
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



## ## ## ## ## ## ## ## ADD BATTERIES ## ## ## ## ## ## ## ##
# annualizing capital (overnight and FOM) costs for battery (inverter) power 
# capacity over 20 years with a discount rate of 0.07
capital_cost_inverter = annuity(20,0.07)*310e3*(1 + 0.01) # €/MW

# capital cost of batteries for every MWh storage capacity, annualized over 
# 15 years with a discount rate of 0.07
capital_cost_batteries = annuity(20,0.07)*144.6e3*(1 + 0) # €/MWh

# battery inverter efficiency
efficiency_inverter = 0.9


# adding the Austrian battery bus
network.add("Bus",
            "Austrian battery bus",
            carrier="DC")

# adding the Austrian link to charge and discharge the batteries from the grid
network.add("Link",
            "Austrian inverter",
            bus0="Austrian battery bus",
            bus1="Austrian electricity bus",
            p_nom_extendable=True,
            p_min_pu=-1, # link is reversible
            capital_cost=capital_cost_inverter*efficiency_inverter,
            marginal_cost=0,
            efficiency=efficiency_inverter
            )

# adding the Austrian store to store the electricity in the batteries
network.add("Store",
            "Austrian batteries",
            bus="Austrian battery bus",
            e_nom_extendable=True,
            capital_cost=capital_cost_batteries,
            e_cyclic=True,
            )


# adding the Slovak battery bus
network.add("Bus",
            "Slovak battery bus",
            carrier="DC")

# adding the Slovak link to charge and discharge the batteries from the grid
network.add("Link",
            "Slovak inverter",
            bus0="Slovak battery bus",
            bus1="Slovak electricity bus",
            p_nom_extendable=True,
            p_min_pu=-1, # link is reversible
            capital_cost=capital_cost_inverter*efficiency_inverter,
            marginal_cost=0,
            efficiency=efficiency_inverter
            )

# adding the Slovak store to store the electricity in the batteries
network.add("Store",
            "Slovak batteries",
            bus="Slovak battery bus",
            e_nom_extendable=True,
            capital_cost=capital_cost_batteries,
            e_cyclic=True,
            )


# adding the Czech battery bus
network.add("Bus",
            "Czech battery bus",
            carrier="DC")

# adding the Czech link to charge and discharge the batteries from the grid
network.add("Link",
            "Czech inverter",
            bus0="Czech battery bus",
            bus1="Czech electricity bus",
            p_nom_extendable=True,
            p_min_pu=-1, # link is reversible
            capital_cost=capital_cost_inverter*efficiency_inverter,
            marginal_cost=0,
            efficiency=efficiency_inverter
            )

# adding the Czech store to store the electricity in the batteries
network.add("Store",
            "Czech batteries",
            bus="Czech battery bus",
            e_nom_extendable=True,
            capital_cost=capital_cost_batteries,
            e_cyclic=True,
            )



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


# adding the Austrian hydrogen bus
network.add("Bus",
            "Austrian hydrogen bus",
            carrier="hydrogen")

# adding the Austrian link to act as the fuel cell
network.add("Link",
            "Austrian fuel cell",
            bus0="Austrian hydrogen bus",
            bus1="Austrian electricity bus",
            p_nom_extendable=True,
            capital_cost=capital_cost_fuelCell*efficiency_fuelCell,
            efficiency=efficiency_fuelCell
            )

# adding the Austrian link to act as the elctrolysis
network.add("Link",
            "Austrian electrolysis",
            bus0="Austrian electricity bus",
            bus1="Austrian hydrogen bus",
            p_nom_extendable=True,
            capital_cost=capital_cost_electrolysis*efficiency_electrolysis,
            efficiency=efficiency_electrolysis
            )

# adding the Austrian store to store the hydrogen
network.add("Store",
            "Austrian hydrogen storage",
            bus="Austrian hydrogen bus",
            e_nom_extendable=True,
            capital_cost=capital_cost_hydrogen,
            e_cyclic=True
            )


# adding the Slovak hydrogen bus
network.add("Bus",
            "Slovak hydrogen bus",
            carrier="hydrogen")

# adding the Slovak link to act as the fuel cell
network.add("Link",
            "Slovak fuel cell",
            bus0="Slovak hydrogen bus",
            bus1="Slovak electricity bus",
            p_nom_extendable=True,
            capital_cost=capital_cost_fuelCell*efficiency_fuelCell,
            efficiency=efficiency_fuelCell
            )

# adding the Slovak link to act as the elctrolysis
network.add("Link",
            "Slovak electrolysis",
            bus0="Slovak electricity bus",
            bus1="Slovak hydrogen bus",
            p_nom_extendable=True,
            capital_cost=capital_cost_electrolysis*efficiency_electrolysis,
            efficiency=efficiency_electrolysis
            )

# adding the Slovak store to store the hydrogen
network.add("Store",
            "Slovak hydrogen storage",
            bus="Slovak hydrogen bus",
            e_nom_extendable=True,
            capital_cost=capital_cost_hydrogen,
            e_cyclic=True
            )


# adding the Czech hydrogen bus
network.add("Bus",
            "Czech hydrogen bus",
            carrier="hydrogen")

# adding the Czech link to act as the fuel cell
network.add("Link",
            "Czech fuel cell",
            bus0="Czech hydrogen bus",
            bus1="Czech electricity bus",
            p_nom_extendable=True,
            capital_cost=capital_cost_fuelCell*efficiency_fuelCell,
            efficiency=efficiency_fuelCell
            )

# adding the Czech link to act as the elctrolysis
network.add("Link",
            "Czech electrolysis",
            bus0="Czech electricity bus",
            bus1="Czech hydrogen bus",
            p_nom_extendable=True,
            capital_cost=capital_cost_electrolysis*efficiency_electrolysis,
            efficiency=efficiency_electrolysis
            )

# adding the Czech store to store the hydrogen
network.add("Store",
            "Czech hydrogen storage",
            bus="Czech hydrogen bus",
            e_nom_extendable=True,
            capital_cost=capital_cost_hydrogen,
            e_cyclic=True
            )



## ## ## ## ## ## ## ## CONNECT THE NETWORKS ## ## ## ## ## ## ## ##
# distance between capitals [km]
dist_AUT_SVK = 51.54
dist_AUT_CZE = 256.31
dist_CZE_SVK = 289.00

# overnight costs for HVDC lines and a converter pair [€/MW]
ON_AUT_SVK = 400*dist_AUT_SVK + 150e3
ON_AUT_CZE = 400*dist_AUT_CZE + 150e3
ON_CZE_SVK = 400*dist_CZE_SVK + 150e3

# capital cost for HVDC lines + converter pair annualized over 40 years with a 
# discount rate of 0.07 [€/MW]
capital_cost_AUT_SVK = annuity(40,0.07)*ON_AUT_SVK*(1 + 0.02)
capital_cost_AUT_CZE = annuity(40,0.07)*ON_AUT_CZE*(1 + 0.02)
capital_cost_CZE_SVK = annuity(40,0.07)*ON_CZE_SVK*(1 + 0.02)

# adding links between the capitals
network.add("Link",
            "AUT - SVK",
            bus0="Austrian electricity bus",
            bus1="Slovak electricity bus",
            p_nom_extendable=True, # transmission capacity is optimized
            p_min_pu=-1, # the link is reversible
            length=dist_AUT_SVK,
            capital_cost=capital_cost_AUT_SVK
            )
network.add("Link",
            "AUT - CZE",
            bus0="Austrian electricity bus",
            bus1="Czech electricity bus",
            p_nom_extendable=True, # transmission capacity is optimized
            p_min_pu=-1, # the link is reversible
            length=dist_AUT_CZE,
            capital_cost=capital_cost_AUT_CZE
            )
network.add("Link",
            "CZE - SVK",
            bus0="Czech electricity bus",
            bus1="Slovak electricity bus",
            p_nom_extendable=True, # transmission capacity is optimized
            p_min_pu=-1, # the link is reversible
            length=dist_CZE_SVK,
            capital_cost=capital_cost_CZE_SVK
            )



## ## ## ## ## ## ## ## SOLVE THE NETWORK ## ## ## ## ## ## ## ##
network.lopf(network.snapshots,solver_name="gurobi")


#%%
## ## ## ## ## ## ## ## DATA-PROCESSING ## ## ## ## ## ## ## ##
# plotting the annual electricity mix of Austria
fig = plt.figure()
labels = ["hydro","onshore wind", "solar", "gas (OCGT)"]
sizes = [network.storage_units_t.p["Austrian hydro"].sum(),
         network.generators_t.p["Austrian onshorewind"].sum(),
         network.generators_t.p["Austrian solar"].sum(),
         network.generators_t.p["Austrian OCGT"].sum()]
colors=["blue", "black", "orange", "brown"]

plt.pie(sizes,
        colors=colors,
        labels=labels,
        wedgeprops={"linewidth":0})
plt.axis("equal")
plt.title("Austria")
plt.tight_layout()
plt.savefig("../LaTeX/figures/G_mix_AUT.eps")

# plotting the annual electricity mix of Slovakia
fig = plt.figure()
labels = ["hydro","onshore wind", "solar", "gas (OCGT)"]
sizes = [network.storage_units_t.p["Slovak hydro"].sum(),
         network.generators_t.p["Slovak onshorewind"].sum(),
         network.generators_t.p["Slovak solar"].sum(),
         network.generators_t.p["Slovak OCGT"].sum()]
colors=["blue", "black", "orange", "brown"]

plt.pie(sizes,
        colors=colors,
        labels=labels,
        wedgeprops={"linewidth":0})
plt.axis("equal")
plt.title("Slovakia")
plt.tight_layout()
plt.savefig("../LaTeX/figures/G_mix_SVK.eps")

# plotting the annual electricity mix of Czechia
fig = plt.figure()
labels = ["hydro","onshore wind", "solar", "gas (OCGT)"]
sizes = [network.storage_units_t.p["Czech hydro"].sum(),
         network.generators_t.p["Czech onshorewind"].sum(),
         network.generators_t.p["Czech solar"].sum(),
         network.generators_t.p["Czech OCGT"].sum()]
colors=["blue", "black", "orange", "brown"]

plt.pie(sizes,
        colors=colors,
        labels=labels,
        wedgeprops={"linewidth":0})
plt.axis("equal")
plt.title("Czechia")
plt.tight_layout()
plt.savefig("../LaTeX/figures/G_mix_CZE.eps")


# timeperiod of dispatch time series, three days in summer
startdate = "2015-06-07"
enddate = "2015-06-10"

# plotting the dispatch time series for Austria
fig = plt.figure()
plt.plot(network.loads_t.p["Austrian load"].loc[startdate:enddate], 
         color="black", linestyle="dashed", label="demand")
plt.plot(network.storage_units_t.p["Austrian hydro"].loc[startdate:enddate], 
         color="blue", label="hydro")
plt.plot(network.generators_t.p["Austrian onshorewind"].loc[startdate:enddate], 
         color="black", label="onshore wind")
plt.plot(network.generators_t.p["Austrian solar"].loc[startdate:enddate], 
         color="orange", label="solar")
plt.plot(network.generators_t.p["Austrian OCGT"].loc[startdate:enddate], 
         color="brown", label="gas (OCGT)")
plt.plot(network.links_t.p1["Austrian inverter"].loc[startdate:enddate], 
         color="cyan", label="batteries")
plt.plot((-network.links_t.p1["Austrian fuel cell"] 
          -network.links_t.p0["Austrian electrolysis"]).loc[startdate:enddate], 
         color="purple", label="hydrogen")
fig.autofmt_xdate()
plt.xlim([pd.Timestamp(startdate), pd.Timestamp(enddate)])
plt.ylim([-15000,30000])
plt.ylabel("MW electricity")
plt.grid(True)
plt.legend(loc="center", bbox_to_anchor=(.5,1.1), frameon=False , ncol=4)
plt.tight_layout()
plt.savefig("../LaTeX/figures/G_disp_AUT.eps")

# plotting the dispatch time series for Slovakia
fig = plt.figure()
plt.plot(network.loads_t.p["Slovak load"].loc[startdate:enddate], 
         color="black", linestyle="dashed", label="demand")
plt.plot(network.storage_units_t.p["Slovak hydro"].loc[startdate:enddate], 
         color="blue", label="hydro")
plt.plot(network.generators_t.p["Slovak onshorewind"].loc[startdate:enddate], 
         color="black", label="onshore wind")
plt.plot(network.generators_t.p["Slovak solar"].loc[startdate:enddate], 
         color="orange", label="solar")
plt.plot(network.generators_t.p["Slovak OCGT"].loc[startdate:enddate], 
         color="brown", label="gas (OCGT)")
plt.plot(network.links_t.p1["Slovak inverter"].loc[startdate:enddate], 
         color="cyan", label="batteries")
plt.plot((-network.links_t.p1["Slovak fuel cell"] 
          -network.links_t.p0["Slovak electrolysis"]).loc[startdate:enddate], 
         color="purple", label="hydrogen")
fig.autofmt_xdate()
plt.xlim([pd.Timestamp(startdate), pd.Timestamp(enddate)])
plt.ylim([-1000, 4000])
plt.ylabel("MW electricity")
plt.grid(True)
plt.legend(loc="center", bbox_to_anchor=(.5,1.1), frameon=False , ncol=4)
plt.tight_layout()
plt.savefig("../LaTeX/figures/G_disp_SVK.eps")

# plotting the dispatch time series for Czechia
fig = plt.figure()
plt.plot(network.loads_t.p["Czech load"].loc[startdate:enddate], 
         color="black", linestyle="dashed", label="demand")
plt.plot(network.storage_units_t.p["Czech hydro"].loc[startdate:enddate], 
         color="blue", label="hydro")
plt.plot(network.generators_t.p["Czech onshorewind"].loc[startdate:enddate], 
         color="black", label="onshore wind")
plt.plot(network.generators_t.p["Czech solar"].loc[startdate:enddate], 
         color="orange", label="solar")
plt.plot(network.generators_t.p["Czech OCGT"].loc[startdate:enddate], 
         color="brown", label="gas (OCGT)")
plt.plot(network.links_t.p1["Czech inverter"].loc[startdate:enddate], 
         color="cyan", label="batteries")
plt.plot((-network.links_t.p1["Czech fuel cell"] 
          -network.links_t.p0["Czech electrolysis"]).loc[startdate:enddate], 
         color="purple", label="hydrogen")
fig.autofmt_xdate()
plt.xlim([pd.Timestamp(startdate), pd.Timestamp(enddate)])
plt.ylim([-2000, 10000])
plt.ylabel("MW electricity")
plt.grid(True)
plt.legend(loc="center", bbox_to_anchor=(.5,1.1), frameon=False , ncol=4)
plt.tight_layout()
plt.savefig("../LaTeX/figures/G_disp_CZE.eps")


print("\nNominal wind power capacity in Austria: {:4.0f} MW".format(
    network.generators.p_nom_opt["Austrian onshorewind"]))
print("Nominal solar power capacity in Austria: {:4.0f} MW".format(
    network.generators.p_nom_opt["Austrian solar"]))
print("Nominal OCGT power capacity in Austria: {:4.0f} MW".format(
    network.generators.p_nom_opt["Austrian OCGT"]))
print("Nominal hydro power capacity in Austria: {:4.0f} MW".format(
    network.storage_units.p_nom_opt["Austrian hydro"]))

print("\nNominal battery power capacity in Austria: {:4.0f} MW".format(
    network.links.p_nom_opt["Austrian inverter"]))
print("Nominal battery storage capacity in Austria: {:4.0f} MWh".format(
    network.stores.e_nom_opt["Austrian batteries"]))

print("\nNominal electrolysis power capacity in Austria: {:4.0f} MW".format(
    network.links.p_nom_opt["Austrian electrolysis"]))
print("Nominal fuel cell power capacity in Austria: {:4.0f} MW".format(
    network.links.p_nom_opt["Austrian fuel cell"]))
print("Nominal hydrogen storage capacity in Austria: {:4.0f} MWh".format(
    network.stores.e_nom_opt["Austrian hydrogen storage"]))


print("\n\nNominal wind power capacity in Slovakia: {:4.0f} MW".format(
    network.generators.p_nom_opt["Slovak onshorewind"]))
print("Nominal solar power capacity in Slovakia: {:4.0f} MW".format(
    network.generators.p_nom_opt["Slovak solar"]))
print("Nominal OCGT power capacity in Slovakia: {:4.0f} MW".format(
    network.generators.p_nom_opt["Slovak OCGT"]))
print("Nominal hydro power capacity in Slovakia: {:4.0f} MW".format(
    network.storage_units.p_nom_opt["Slovak hydro"]))

print("\nNominal battery power capacity in Slovakia: {:4.0f} MW".format(
    network.links.p_nom_opt["Slovak inverter"]))
print("Nominal battery storage capacity in Slovakia: {:4.0f} MWh".format(
    network.stores.e_nom_opt["Slovak batteries"]))

print("\nNominal electrolysis power capacity in Slovakia: {:4.0f} MW".format(
    network.links.p_nom_opt["Slovak electrolysis"]))
print("Nominal fuel cell power capacity in Slovakia: {:4.0f} MW".format(
    network.links.p_nom_opt["Slovak fuel cell"]))
print("Nominal hydrogen storage capacity in Slovakia: {:4.0f} MWh".format(
    network.stores.e_nom_opt["Slovak hydrogen storage"]))


print("\n\nNominal wind power capacity in Czechia: {:4.0f} MW".format(
    network.generators.p_nom_opt["Czech onshorewind"]))
print("Nominal solar power capacity in Czechia: {:4.0f} MW".format(
    network.generators.p_nom_opt["Czech solar"]))
print("Nominal OCGT power capacity in Czechia: {:4.0f} MW".format(
    network.generators.p_nom_opt["Czech OCGT"]))
print("Nominal hydro power capacity in Czechia: {:4.0f} MW".format(
    network.storage_units.p_nom_opt["Czech hydro"]))

print("\nNominal battery power capacity in Czechia: {:4.0f} MW".format(
    network.links.p_nom_opt["Czech inverter"]))
print("Nominal battery storage capacity in Czechia: {:4.0f} MWh".format(
    network.stores.e_nom_opt["Czech batteries"]))

print("\nNominal electrolysis power capacity in Czechia: {:4.0f} MW".format(
    network.links.p_nom_opt["Czech electrolysis"]))
print("Nominal fuel cell power capacity in Czechia: {:4.0f} MW".format(
    network.links.p_nom_opt["Czech fuel cell"]))
print("Nominal hydrogen storage capacity in Czechia: {:4.0f} MWh".format(
    network.stores.e_nom_opt["Czech hydrogen storage"]))


print("\n\nPower capacity of Austrian-Slovak link: {:4.0f} MW".format(
    network.links.p_nom_opt["AUT - SVK"]))
print("Power capacity of Austrian-Czech link:  {:4.0f} MW".format(
    network.links.p_nom_opt["AUT - CZE"]))
print("Power capacity of Czech-Slovak link:     {:3.0f} MW".format(
    network.links.p_nom_opt["CZE - SVK"]))

print("\nCO2-limit: {:4.0f} kt".format(
    network.global_constraints.constant["co2_limit"]*1e-3))
print("CO2-price: {:2.0f} €".format(
    network.global_constraints.mu["co2_limit"]))

print("\nTotal system cost: {:4.0f} million €".format(
    network.objective*1e-6))
print("Marginal system cost: {:2.1f} €/MWh".format( \
    network.objective/(network.loads_t.p["Austrian load"].sum()
                       + network.loads_t.p["Slovak load"].sum()
                       + network.loads_t.p["Czech load"].sum())))
