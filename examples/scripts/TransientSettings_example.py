from corems.encapsulation.factory.parameters import MSParameters
from corems.encapsulation.factory.processingSetting import TransientSetting

print(TransientSetting.apodization_method) # Hanning, as expected bc it is the default value
TransientSetting.apodization_method = 'Full-Sine'
print(TransientSetting.apodization_method) # Full-Sine, as expected bc it was changed on the class

# Making an instance of TransientSetting uses the default value for apodization_method, not the value set on the class above
transient_instance = TransientSetting()
print(transient_instance.apodization_method) # Hanning, used the default value - is this expected?
transient_instance.apodization_method = 'Half-Sine'
print(transient_instance.apodization_method) # Half-Sine, as expected bc it was changed on the instance
print(TransientSetting.apodization_method) # Full-Sine, as expected bc it was changed on the class (but not used by any instance)

print(MSParameters.transient.apodization_method) # Hanning, used the default value, even though it was changed on the TransientSetting class