def normalize_coordinates(coordinates):
    return [coord / max(coordinates) for coord in coordinates]
