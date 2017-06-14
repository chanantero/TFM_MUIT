classdef (ConstructOnLoad) updatedValuesEvntData < event.EventData
   properties
      type
   end
   
   methods
      function data = updatedValuesEvntData(type)
         data.type = type;
      end
   end
end