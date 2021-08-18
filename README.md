# UK-distribution-hub-finder

Hillclimbing-based search algorithm for finding optimal locations for distribution hubs in the UK.
Based on the locations and populations of top 100 cities in UK, this program implements an optimised hill-climb search for the optimal location for a national distribution hub, as well as its nearest neighbouring city.

There are a number of customisable parameters, including:
 - including the salaries of the nearest cities, which would impact the operating costs of the distribution centre
 - including the location of ports in the UK as a factor influencing the position of a distribution hub
 - choosing the number of distribution hubs for which their locations need optimising (decided by a custom _k_-means clustering algorithm)

The program returns the location of all hubs, the nearest cities to these hubs, the runtime of the program and the number of calculations needed during optimisation.

Note, that the distance between locations is assumed to be as the crow flies, using the [Haversine formula](https://en.wikipedia.org/wiki/Haversine_formula) for distances on Earth.
As such, road mileage between cities and time of travel are not currently implemented for optimising the position of hubs.
