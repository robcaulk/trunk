#include "BodyContainer.hpp"
#include "Body.hpp"

void BodyContainer::registerAttributes()
{
	REGISTER_ATTRIBUTE(body);
};

void BodyContainer::preProcessAttributes(bool deserializing)
{
	if(deserializing)
	{
		body.clear();
	}
	else
	{
		body.clear();
		for( this->gotoFirst() ; this->notAtEnd() ; this->gotoNext() )
			body.push_back(this->getCurrent());
	}
};

void BodyContainer::postProcessAttributes(bool deserializing)
{
	if(deserializing)
	{
		this->clear();
		vector<shared_ptr<Body> >::iterator it    = body.begin();
		vector<shared_ptr<Body> >::iterator itEnd = body.end();
		for( ; it != itEnd ; ++it)
			this->insert(*it);
		body.clear();
	}
	else
	{
		body.clear();
	}
};

void BodyContainer::setId(shared_ptr<Body>& b, unsigned int newId)
{
	b->id = newId;
}
