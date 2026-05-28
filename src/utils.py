







def parse_poni(poni_str):
    raw_data = {}
    for line in poni_str.strip().split("\n"):
       
        if line.startswith("#") or not line:
            continue
        key, value = line.split(":", 1)
        raw_data[key.strip()] = value.strip()
        
    return  float(raw_data["Poni1"]),float(raw_data["Poni2"]),
    
    