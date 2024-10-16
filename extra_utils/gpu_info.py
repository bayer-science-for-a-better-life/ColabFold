import jax
import jaxlib
import subprocess
import tensorflow as tf

# Check CUDA version if available
try:
    from jax.lib import xla_bridge

    print("CUDA version:", xla_bridge.get_backend().platform)
except ImportError:
    print("CUDA is not available")

# Check cuDNN version if available
try:
    from jax.lib import cudnn

    print("cuDNN version:", cudnn.getVersion())
except ImportError:
    print("cuDNN is not available")

# Check CUDA Toolkit version if available
try:
    from jax.lib import cuda

    print("CUDA Toolkit version:", cuda.cudart.getVersion())
except ImportError:
    print("CUDA Toolkit is not available")


def get_nvcc_version():
    try:
        result = subprocess.run(['nvcc', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"nvcc command failed with error: {result.stderr}")

        # Extract version information from the output
        output_lines = result.stdout.split('\n')
        for line in output_lines:
            if 'release' in line:
                version_info = line.strip()
                return version_info
        raise RuntimeError("Could not find version information in nvcc output")
    except FileNotFoundError:
        return "nvcc is not installed or not found in PATH"


# Example usage
if __name__ == "__main__":
    print("#" * 80)
    print("TensorFlow version:", tf.__version__)
    print("JAX version:", jax.__version__)
    print("JAXlib version:", jaxlib.__version__)
    print(f'JAX backend: {jax.default_backend()}')
    nvcc_version = get_nvcc_version()
    print("nvcc version:", nvcc_version)
