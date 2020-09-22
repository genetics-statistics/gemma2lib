from types import SimpleNamespace

def methodize(data: dict) -> SimpleNamespace:
    """Make items callable as a method - this is not the final version
    because we want to make it read-only

    """
    return SimpleNamespace(**data)
