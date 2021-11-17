#ifndef DGG_MODEL_HPP
#define DGG_MODEL_HPP

namespace Cajete 
{

//Simple specification for a DGG Model base type
//
//It's template on a generic interface type, this 
//is inteded to ensure the model can be set by any
//command line parser
template <typename InterfaceType>
class DggModel {
    public:
        using interface_type = InterfaceType;

        virtual ~DggModel() = default;

        virtual void init(InterfaceType interface) const = 0;

        virtual void run() const = 0;
};

} //end namespace Cajete

#endif 
