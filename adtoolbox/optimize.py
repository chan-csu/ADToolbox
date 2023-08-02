import openbox
import adm

class Tuner:
    """
    This class is a wrapper around openbox to optimize a model.
    """
    def __init__(self,
                 base_model: adm.Model,
                 train_data: list[dict],
                 fitness_mode:str="equalized",
                 **kwargs)->None:
        """
        base_model: The model to optimize.
        train_data: The data to train the model on.
        fitness_mode: The fitness mode to use.
        kwargs: Additional arguments to pass to openbox.
        """
        self.base_model = base_model
        self.train_data = train_data
        self.fitness_mode = fitness_mode
        self.kwargs = kwargs
        
    def fitness(self, config: dict)->float:
        """
        This function is called by openbox to evaluate a configuration.
        :param config: The configuration to evaluate.
        :return: The fitness of the configuration.
        """
        if self.fitness_mode == "equalized":
            for experiment in self.train_data:
                self._model= self.base_model.copy()
                
            