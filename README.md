# ❀ Swarm Intelligence in Bioinformatics ❀

### Particle Swarm Optimization (PSO)


### Artificial Bee Colony (ABC)

### Ant Colony Optimization (ACO)
Artificial Bee Colony was first proposed by Karaboga in 2005. It is an optimization algorithm inspired by the intelligent foraging behavior of a honey bee swarm for finding an optimal solution. Three essential components of forage selection are:
1)	food sources: places with a high amount of nectar 
2)	employed foragers: employed bees exploit the food source, calculate the nector amount, and carry the optimal path information back to hive.
3)	unemployed foragers: unemployed bees consist of two groups of bees – onlookers and scouts. Onlookers wait in the hive for the information that is shared by the employed bees. Scout bees are translated from employed bees and always continues seeking new food sources near the hive.   
Comparing Artificial Bee Colony algorithm with Particle Swarm Optimization and Ant Colony Optimization methods, it carries the advantages of fast convergence, high flexibility, and strong robustness. The mechanism of ABC is provided in Algorithm 3.
```
Algorithm 3: Artificial Bee Colony
1:	Initialize food sources
2:	while termination criteria is not met do
3:	for each employed bee
4:	Produce new solution
5:	Calculate the fitness value
6:	Apply greedy selection 
7:	Calculate the probability value
8:	end
9:	for each onlooker bee
10:	Select a solution 
11:	Produce new solution
12:	Calculate the fitness value
13:	Apply greedy selection 
14:	end
15:	if an abandoned solution for the scout exists then
16:	Replace it with a new solution at random
17:	end
18:	Register the best solution
19:	end
```

