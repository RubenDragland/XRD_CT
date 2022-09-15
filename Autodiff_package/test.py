
import torch
x = torch.tensor([0.5, 0.75], requires_grad=True)
y = torch.log(x[0] * x[1]) * torch.sin(x[1])
y.backward()
print(x.grad)

