CellType = c("OL", "OL", "OL", "OL", "OL", "OL", "OL",
             "Neuron", "Neuron", "Astro", "Neuron", "Astro", "Neuron", "Astro",
             "Neuron", "Neuron", "Astro", "Neuron", "Astro", "Astro", "Astro",
             "Astro", "Astro", "OL", "OL", "OL", "Astro", "Astro", "Neuron",
             "Astro", "Astro", "Astro", "Astro")

Contrasts = makeContrasts(NeuronOL = CellTypeNeuron - CellTypeOL,
                          NeuronAstro = CellTypeNeuron - CellTypeAstro,
                          OLAstro = CellTypeOL - CellTypeAstro,
                          levels = Design)