import re
import os

"""
Written with significant generation from Claude Sonnet 4.5
"""



def get_n_used(date, target_name, base_dir):
    """
    Count the number of lines in mag<date>_<target_name>.lis
    
    Args:
        date: Date as a string
        target_name: Target name as a string
        base_dir: Base directory path
    
    Returns:
        N_used: Number of lines in the file
    """
    filename = os.path.join(base_dir + '/combo/', f"mag{date}_{target_name}.lis")
    
    try:
        with open(filename, 'r') as f:
            n_used = sum(1 for line in f)
        return n_used
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return None


def get_ref_and_error(target_name, base_dir):
    """
    Extract N_ref and Ast_Error from plotPosError_<target>.txt
    
    Args:
        target_name: Target name as a string
        base_dir: Base directory path
    
    Returns:
        tuple: (N_ref, Ast_Error) or (None, None) if error
    """
    filename = os.path.join(base_dir + '/combo/starfinder/', f"plotPosError_{target_name}.txt")
    
    try:
        with open(filename, 'r') as f:
            content = f.read()
        
        # Extract number of detections
        n_ref_match = re.search(r'Number of detections:\s+(\d+)', content)
        n_ref = int(n_ref_match.group(1)) if n_ref_match else None
        
        # Extract median position error
        ast_error_match = re.search(r'Median Pos Error \(mas\).*?:\s+([\d.]+)', content)
        ast_error = float(ast_error_match.group(1)) if ast_error_match else None
        
        return n_ref, ast_error
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return None, None

def get_strehl_and_fwhm(date, target_name, base_dir):
    """
    Extract Strehl and FWHM from stf_psf_strehl_<target>.txt
    
    Args:
        date: Date as a string
        target_name: Target name as a string
        base_dir: Base directory path
    
    Returns:
        tuple: (strehl, fwhm) or (None, None) if error
    """
    filename = os.path.join(base_dir + '/combo/starfinder/', f"stf_psf_strehl_{target_name}.txt")
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                # Look for the line with mag<date>_<target>_psf.fits
                if f"mag{date}_{target_name}_psf.fits" in line:
                    # Split the line by whitespace
                    parts = line.split()
                    if len(parts) >= 4:
                        strehl = float(parts[1])
                        fwhm = float(parts[3])
                        return strehl, fwhm
        
        print(f"Error: Could not find matching line in '{filename}'.")
        return None, None
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return None, None
    except (ValueError, IndexError) as e:
        print(f"Error parsing values: {e}")
        return None, None

def get_astrometric_errors(date, target_name, object_name, base_dir):
    """
    Extract x_ast_err and y_ast_err from mag<date>_<target>_rms.lis
    
    Args:
        date: Date as a string
        target_name: Target name as a string
        object_name: Name of the object to search for in the file
        base_dir: Base directory path
    
    Returns:
        tuple: (x_ast_err, y_ast_err) or (None, None) if error
    """
    filename = os.path.join(base_dir + '/combo/starfinder/', f"mag{date}_{target_name}_rms.lis")
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                # Skip comment lines
                if line.startswith('#'):
                    continue
                
                # Split the line by whitespace
                parts = line.split()
                
                # Check if this is the line for the specified object
                if len(parts) > 0 and parts[0] == object_name:
                    if len(parts) >= 6:
                        x_ast_err = float(parts[5])  # xe is the 5th column (index 4)
                        y_ast_err = float(parts[6])  # ye is the 6th column (index 5)
                        return x_ast_err, y_ast_err
        
        print(f"Error: Could not find object '{object_name}' in '{filename}'.")
        return None, None
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return None, None
    except (ValueError, IndexError) as e:
        print(f"Error parsing values: {e}")
        return None, None

def get_summary_stats(date, target_name, object_name, base_dir):
    """
    date : str
       Epoch name such as 21jun13os
    target_name : str
       Target name with suffix such as m15_kp
    object_name : str
       Object name as it appears in starfinder, such as m15
    base_dir : str
       Base directory such as /g/lu/data/microlens/21jun13os/
    """

    n_used = get_n_used(date, target_name, base_dir)
    n_used = get_n_used(date, target_name, base_dir)
    n_ref, ast_error = get_ref_and_error(target_name, base_dir)
    strehl, fwhm = get_strehl_and_fwhm(date, target_name, base_dir)
    x_ast_err, y_ast_err = get_astrometric_errors(date, target_name, object_name, base_dir)

    if all(v is not None for v in [n_used, n_ref, ast_error, strehl, fwhm, x_ast_err, y_ast_err]):
        # Print header row
        print("N_used\tstrehl\tfwhm\tAst_Error\tN_ref\tx_ast_err\ty_ast_err")
        # Print values row
        print(f"{n_used}\t{strehl}\t{fwhm}\t{ast_error}\t{n_ref}\t{x_ast_err}\t{y_ast_err}")

    return
