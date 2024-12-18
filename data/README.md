# Data for replication
This file defines the various datasets used in the paper, specifically regarding the arrival rate of the customers, and details on the various scheduling policies used in the paper.

## Arrival Rates
The ArrivalRates.xlsx file includes the mean arrival rates (in arrivals/hr) for first-time and crisis services across the 18 weeks of the semester. 

## Schedule proportions under various policies obtained by solving **S-WST** 
The ScheduleProportions.xlsx file includes the schedule proportions obtained by solving Problem **S-WST** for first-time and crisis services across the 18 weeks of the semester, obtained using the methodology proposed in the. paper It includes 18 sheets, with each sheet name in the format {F}-{C}, where {F} represents the average overall first-time service proportion, and {C} represents the average overall crisis service proportion in the schedule. For instance, the sheet "8-4" represents the case when the first-time service proportion is set to 8% and the crisis service proportion is set to 4%. In the sheet "8-4", the value of 0.125 in Week 1 for first-time service indicates that across the week, 12.5% of the total number of slots were dedicated to serving first-time customers, on average. These proportions are used to generate the results in Table 1 and Figure 9 of the paper.

## Schedule proportions under the DPS policy
The DPS_Policy.xlsx file includes the scheduling proportions obtained using the DPS policy, shown in Figure 10 and used to generate the results of Figure 11 of the paper.
