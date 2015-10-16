/* Siconos-Kernel, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

/*! \file EventFactory.hpp
\brief  Factory to generate user-defined Events
*/

#ifndef EventFactory_H
#define EventFactory_H

# include "Event.hpp"
# include <map>


/** Namespace for Events factory related objects. */
namespace EventFactory
{

/** A pointer to function, returning a pointer to Event, built with its type (ie class name) and a pointer to Model.*/
typedef SP::Event(*object_creator)(double, int);

/** The type of the factory map */
typedef std::map<int, object_creator> MapFactory;

/** An iterator through the MapFactory */
typedef MapFactory::iterator MapFactoryIt;

/** Template function to return a new object of type SubType
 * \param time time of the Event
 * \param type type of the Event
 * \return an Event
 */
template<class SubType> SP::Event factory(double time, int type)
{
  SP::Event e(new SubType(time, type));
  return e;
}

/** Registry Class for Events.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) February 01, 2007
 *
 */
class Registry
{

private :

  /** map that links a std::string, the type of the class, to a pointer to function, used to build the object. */
  MapFactory factory_map;

public :

  /** get access to the Registry
   * \return reference to the registry
   */
  static Registry& get() ;

  /** Add an object_creator into the factory_map, factory_map[name] = object.
   * \param type the type of the object added
   * \param creator object creator
   */
  void add(int type, object_creator creator);

  /**
   *  \param time time of Event
   *  \param type type of Event
   *  \return an Event
   */
  SP::Event instantiate(double time, int type);

} ;

/** Registration Class for Events.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) February 01, 2007
 *
 * Class used for auto-registration of Event-type objects.
 *
 */
class Registration
{

public :

  /** To register some new object into the factory
   * \param type the type of the object added
   * \param creator object creator
   */
  Registration(int type, object_creator creator);
} ;

}
// end of namespace EventFactory

#define AUTO_REGISTER_EVENT(class_name,class_type) EventFactory::Registration _registration_## class_type(class_name,&EventFactory::factory<class_type>);
#endif












