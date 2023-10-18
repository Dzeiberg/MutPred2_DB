class Mechanism:
    def __init__(self, name, effect_type, p_value, posterior, position=None) -> None:
        self.name = name
        self.effect_type = effect_type
        self.p_value = p_value
        self.posterior = posterior
        self.position = int(position) if position is not None else None
