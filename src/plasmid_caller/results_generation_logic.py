## hit selection logic!

# Function to determine the best match
def get_best_match(matches, key):
    """Highest percent identity takes the cake. specify which feature to compare.
    ex: get_best_match(matches, "percent_identity")"""
    best_match = None
    best_score = -1
    for match in matches:
        if match[key] > best_score:
            best_score = match[key]
            best_match = match
    return best_match
