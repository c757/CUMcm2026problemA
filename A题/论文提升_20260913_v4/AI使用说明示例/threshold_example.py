def drying_finished(values, threshold=0.15):
    values = tuple(values)
    if not values:
        raise ValueError("empty concentration data")
    return max(values) < threshold


state = [0.14996, 0.12000, 0.05270]
finished = drying_finished(state)

# Round for display only, after the decision.
display_values = [round(c, 4) for c in state]
